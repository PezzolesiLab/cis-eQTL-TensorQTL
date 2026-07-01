#!/usr/bin/env python3
"""Run cis-eQTL mapping with tensorQTL.

The input expression BED is assumed to contain expression values that were
already residualized during preprocessing. Consequently, no additional
covariate matrix is supplied to tensorQTL in this script.

Analyses performed
------------------
1. Gene-level cis-QTL permutation mapping.
2. Nominal cis-QTL mapping, with tensorQTL parquet output by chromosome.
3. Conditional cis-QTL mapping to identify independent signals.
4. Per-chromosome and overall run summaries.

"""

import glob
import gzip
import os
import re
import subprocess
import sys

import pandas as pd
import tensorqtl
import torch
from tensorqtl import cis, genotypeio


# ------------------------------- FUNCTIONS ----------------------------------


def configure_r_environment():
    """Populate R_HOME and LD_LIBRARY_PATH when R is available."""
    if "R_HOME" in os.environ:
        return

    try:
        r_home = subprocess.check_output(["R", "RHOME"], text=True).strip()
        os.environ["R_HOME"] = r_home
        os.environ["LD_LIBRARY_PATH"] = (
            r_home + "/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
        )
    except Exception:
        # R/qvalue is optional because tensorQTL can still report empirical
        # cis permutation p-values when qvalue is unavailable.
        pass


def check_qvalue_availability():
    """Return True when the R qvalue package is available through rpy2."""
    try:
        import rpy2.robjects as ro

        ro.r("suppressMessages(library(qvalue))")
        return True
    except Exception as exc:
        print(
            "Warning: qvalue is not available through R/rpy2:",
            exc,
            file=sys.stderr,
        )
        return False


def normalize_chromosome(chromosome):
    """Convert chromosome labels to chr-prefixed tensorQTL-style names."""
    chromosome = str(chromosome).strip().upper()
    chromosome = re.sub(r"^CHR", "", chromosome)

    if chromosome == "23":
        chromosome = "X"
    elif chromosome == "24":
        chromosome = "Y"
    elif chromosome == "25":
        chromosome = "XY"
    elif chromosome in {"26", "M", "MT"}:
        chromosome = "MT"

    return "chr" + chromosome


def write_tsv_gz_keep_cols(df, out_path, required=("phenotype_id",)):
    """Write a gzipped TSV after confirming required identifier columns."""
    if (
        isinstance(df.index, pd.Index)
        and df.index.name in required
        and df.index.name not in df.columns
    ):
        df = df.reset_index()

    if "phenotype_id" not in df.columns:
        for candidate in ("gene_id", "phenotype", "gene", "Unnamed: 0"):
            if candidate in df.columns:
                df = df.rename(columns={candidate: "phenotype_id"})
                break

    for column in required:
        if column not in df.columns:
            raise KeyError(
                f"Required column '{column}' is missing before writing. "
                f"Present columns: {list(df.columns)}"
            )

    df.to_csv(out_path, sep="\t", index=False, compression="gzip")
    return out_path


def ensure_phenotype_id_column(df):
    """Ensure that an in-memory DataFrame contains a phenotype_id column."""
    if "phenotype_id" not in df.columns:
        if "phenotype" in df.columns:
            df.rename(columns={"phenotype": "phenotype_id"}, inplace=True)
        elif "gene_id" in df.columns:
            df.rename(columns={"gene_id": "phenotype_id"}, inplace=True)

    if df.index.name == "phenotype_id" and "phenotype_id" not in df.columns:
        df.reset_index(inplace=True)

    return df


def chromosome_sort_key(chromosome):
    """Return a natural sorting key for autosomes and sex chromosomes."""
    chromosome = str(chromosome).upper().replace("CHR", "")
    chromosome_order = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}

    try:
        return (0, int(chromosome))
    except ValueError:
        return (1, chromosome_order.get(chromosome, 999))


def bh_fdr(pvalues):
    """Calculate Benjamini-Hochberg adjusted p-values."""
    import numpy as np

    pvalues = np.asarray(pvalues, dtype=float)
    n_tests = pvalues.size
    order = np.argsort(pvalues)
    ranked_pvalues = pvalues[order]
    qvalues = ranked_pvalues * n_tests / (1 + np.arange(n_tests))
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    output = np.empty_like(qvalues)
    output[order.argsort()] = qvalues
    return np.clip(output, 0, 1)


# ------------------------------ CONFIGURATION -------------------------------

plink_prefix_path = "T2T_linear/maf05_allchr_T2Tlinear"
expression_bed = (
    "T2T_norm_residualized_counts_expression.bed.gz"
)
prefix = "T2Tlinear"
outdir = "t2tlinear_tensorqtl_results"

os.makedirs(outdir, exist_ok=True)


# --------------------------- ENVIRONMENT DETAILS ----------------------------

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(
    f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}"
)
print(f"pandas: {pd.__version__}")

configure_r_environment()
have_qvalue = check_qvalue_availability()


# ------------------------------- LOAD INPUTS --------------------------------

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = None
print(
    "Phenotypes loaded without an additional covariate matrix; expression "
    "values are already residualized.",
    file=sys.stderr,
)

plink_reader = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = plink_reader.load_genotypes()
variant_df = plink_reader.bim.set_index("snp")[["chrom", "pos"]]
print("Genotypes loaded.", file=sys.stderr)

phenotype_pos_df = phenotype_pos_df.copy()
phenotype_pos_df.loc[:, "chr"] = phenotype_pos_df["chr"].map(
    normalize_chromosome
)

variant_df = variant_df.copy()
variant_df.loc[:, "chrom"] = variant_df["chrom"].map(normalize_chromosome)

print(
    "Phenotype chromosomes:",
    sorted(phenotype_pos_df["chr"].unique())[:6],
    "...",
)
print(
    "Genotype chromosomes:",
    sorted(variant_df["chrom"].unique())[:6],
    "...",
)


# ---------------------- 1. CIS PERMUTATION MAPPING --------------------------

cis_df = cis.map_cis(
    genotype_df,
    variant_df,
    phenotype_df,
    phenotype_pos_df,
)
ensure_phenotype_id_column(cis_df)

if "phenotype" not in cis_df.columns and "phenotype_id" in cis_df.columns:
    cis_df["phenotype"] = cis_df["phenotype_id"]

try:
    tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
except Exception:
    print(
        "Note: calculate_qvalues() was skipped because R/qvalue was not "
        "available. Empirical permutation p-values are retained.",
        file=sys.stderr,
    )

cis_permutation_file = os.path.join(outdir, f"{prefix}.cis_qtl.txt.gz")
write_tsv_gz_keep_cols(
    cis_df,
    cis_permutation_file,
    required=("phenotype_id",),
)

print(
    "[write] Gene-level cis results:",
    cis_permutation_file,
    file=sys.stderr,
)
print("Cis-QTL permutation mapping complete.", file=sys.stderr)

cis_header = pd.read_csv(cis_permutation_file, sep="\t", nrows=0).columns
assert "phenotype_id" in cis_header


# ------------------------ 2. CIS NOMINAL MAPPING ----------------------------

cis.map_nominal(
    genotype_df,
    variant_df,
    phenotype_df,
    phenotype_pos_df,
    prefix,
    output_dir=outdir,
)
print(
    "[write] Per-chromosome nominal parquet files:",
    outdir,
    file=sys.stderr,
)
print(
    "Nominal cis-QTL mapping completed for all variant-phenotype pairs.",
    file=sys.stderr,
)

parquet_paths = sorted(
    glob.glob(os.path.join(outdir, f"{prefix}.cis_qtl_pairs.chr*.parquet"))
)
nominal_all_file = os.path.join(
    outdir,
    f"{prefix}.cis_qtl_pairs.ALL.tsv.gz",
)
total_nominal_rows = 0

if parquet_paths:
    wrote_header = False
    with gzip.open(nominal_all_file, "wt") as output_handle:
        for file_number, parquet_path in enumerate(parquet_paths, start=1):
            nominal_df = pd.read_parquet(parquet_path, columns=None)
            total_nominal_rows += len(nominal_df)
            nominal_df.to_csv(
                output_handle,
                sep="\t",
                index=False,
                header=not wrote_header,
            )
            wrote_header = True
            print(
                f"[merge] {file_number:02d}/{len(parquet_paths)} -> "
                f"{parquet_path} ({len(nominal_df):,} rows)",
                file=sys.stderr,
            )

    print(
        "[write] Merged nominal TSV:",
        nominal_all_file,
        file=sys.stderr,
    )


# ---------------------- 3. CIS INDEPENDENT SIGNALS --------------------------

if "phenotype_id" in cis_df.columns:
    cis_df = cis_df.set_index("phenotype_id", drop=False)

phenotype_pos_df.index.name = "phenotype_id"

independent_df = cis.map_independent(
    genotype_df,
    variant_df,
    cis_df,
    phenotype_df,
    phenotype_pos_df,
)
ensure_phenotype_id_column(independent_df)

if (
    "phenotype" not in independent_df.columns
    and "phenotype_id" in independent_df.columns
):
    independent_df["phenotype"] = independent_df["phenotype_id"]

independent_file = os.path.join(
    outdir,
    f"{prefix}.cis_independent_qtls.txt.gz",
)
write_tsv_gz_keep_cols(
    independent_df,
    independent_file,
    required=("phenotype_id", "variant_id"),
)

print(
    "[write] Independent cis signals:",
    independent_file,
    independent_df.shape,
    file=sys.stderr,
)
print("Conditional cis-QTL mapping complete.", file=sys.stderr)

independent_header = pd.read_csv(independent_file, sep="\t", nrows=0).columns
assert "phenotype_id" in independent_header
assert "variant_id" in independent_header


# ---------------------- 4. PER-CHROMOSOME SUMMARY ---------------------------

if cis_df.index.name == "phenotype_id":
    cis_df = cis_df.reset_index(drop=True)

if (
    "phenotype_id" in phenotype_pos_df.columns
    and phenotype_pos_df.index.name != "phenotype_id"
):
    phenotype_pos_df = phenotype_pos_df.set_index("phenotype_id", drop=False)

phenotype_pos_df.index.name = "phenotype_id"

variants_by_chromosome = variant_df["chrom"].value_counts().to_dict()
input_phenotypes_by_chromosome = (
    phenotype_pos_df["chr"].value_counts().to_dict()
)

print(
    "DEBUG:",
    "cis_df columns:",
    [column for column in cis_df.columns if column.startswith("phenotype")],
    "| cis_df index name:",
    cis_df.index.name,
    "| phenotype_pos_df index name:",
    phenotype_pos_df.index.name,
    file=sys.stderr,
)

cis_chromosome_df = pd.merge(
    cis_df[["phenotype_id"]].drop_duplicates(),
    phenotype_pos_df[["chr"]],
    left_on="phenotype_id",
    right_index=True,
    how="left",
)
tested_phenotypes_by_chromosome = (
    cis_chromosome_df["chr"].value_counts().to_dict()
)

if "qval" in cis_df.columns:
    egenes_by_chromosome = pd.merge(
        cis_df.loc[cis_df["qval"] <= 0.05, ["phenotype_id"]],
        phenotype_pos_df[["chr"]],
        left_on="phenotype_id",
        right_index=True,
        how="left",
    )["chr"].value_counts().to_dict()
else:
    egenes_by_chromosome = {}

independent_signals_by_chromosome = pd.merge(
    independent_df[["phenotype_id"]].drop_duplicates(),
    phenotype_pos_df[["chr"]],
    left_on="phenotype_id",
    right_index=True,
    how="left",
)["chr"].value_counts().to_dict()

nominal_rows_by_chromosome = {}
for parquet_path in parquet_paths:
    chromosome_match = re.search(
        r"\.chr([A-Za-z0-9]+)\.parquet$",
        parquet_path,
    )
    if not chromosome_match:
        continue

    chromosome = f"chr{chromosome_match.group(1).upper()}"
    try:
        nominal_rows_by_chromosome[chromosome] = len(
            pd.read_parquet(parquet_path, columns=["variant_id"])
        )
    except Exception:
        nominal_rows_by_chromosome[chromosome] = -1

all_chromosomes = (
    set(variants_by_chromosome)
    | set(input_phenotypes_by_chromosome)
    | set(tested_phenotypes_by_chromosome)
    | set(egenes_by_chromosome)
    | set(independent_signals_by_chromosome)
    | set(nominal_rows_by_chromosome)
)

summary_rows = []
for chromosome in sorted(all_chromosomes, key=chromosome_sort_key):
    summary_rows.append(
        {
            "chrom": chromosome,
            "variants": int(variants_by_chromosome.get(chromosome, 0)),
            "phenotypes_input": int(
                input_phenotypes_by_chromosome.get(chromosome, 0)
            ),
            "phenotypes_tested": int(
                tested_phenotypes_by_chromosome.get(chromosome, 0)
            ),
            "egenes_q05": int(egenes_by_chromosome.get(chromosome, 0)),
            "independent_signals": int(
                independent_signals_by_chromosome.get(chromosome, 0)
            ),
            "nominal_pairs_rows": int(
                nominal_rows_by_chromosome.get(chromosome, 0)
            ),
        }
    )

per_chromosome_file = os.path.join(
    outdir,
    f"{prefix}.per_chr_summary.tsv",
)
pd.DataFrame(summary_rows).to_csv(
    per_chromosome_file,
    sep="\t",
    index=False,
)
print(
    "[write] Per-chromosome summary:",
    per_chromosome_file,
    file=sys.stderr,
)


# -------------------------- 5. OVERALL SUMMARY -------------------------------

run_summary = {
    "n_samples": [genotype_df.shape[1]],
    "n_covariates": [0],
    "n_phenotypes_input": [phenotype_df.shape[0]],
    "n_phenotypes_tested": [cis_df.shape[0]],
    "n_egenes_q05": [
        int((cis_df["qval"] <= 0.05).sum()) if "qval" in cis_df.columns else -1
    ],
    "n_independent_signals": [independent_df.shape[0]],
    "parquet_files": [len(parquet_paths)],
    "nominal_pairs_rows_total": [int(total_nominal_rows)],
    "merged_nominal_tsv": [
        os.path.basename(nominal_all_file) if parquet_paths else ""
    ],
    "have_qvalue": [have_qvalue],
}

run_summary_file = os.path.join(outdir, f"{prefix}.run_summary.tsv")
pd.DataFrame(run_summary).to_csv(
    run_summary_file,
    sep="\t",
    index=False,
)
print("[write] Run summary:", run_summary_file, file=sys.stderr)
print("Analysis complete.", file=sys.stderr)
