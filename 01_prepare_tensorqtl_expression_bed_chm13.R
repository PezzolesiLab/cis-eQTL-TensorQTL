#!/usr/bin/env Rscript

# ==============================================================================
# Prepare residualized RNA-seq expression phenotypes for tensorQTL (CHM13/hs1)
# ==============================================================================
#
# Workflow
#   1. Define the final sample set and order from column 2 of a PLINK .fam file.
#   2. Intersect the PLINK samples with RNA-seq count columns and optional keep or
#      exclusion lists.
#   3. Filter genes with edgeR::filterByExpr (or the configured count filter).
#   4. Estimate DESeq2 size factors and apply the variance-stabilizing transform.
#   5. Compute expression principal components and apply a per-gene inverse-
#      normal transformation.
#   6. Residualize expression on configured numeric and categorical covariates.
#   7. Join genes to a precomputed, 0-based CHM13/hs1 TSS table and write a
#      bgzip-compressed, tabix-indexed BED file for tensorQTL.
#
# Expected count-file format
#   The first column contains gene identifiers; all remaining columns are sample
#   IIDs that already match the PLINK IIDs.
#
# Outputs
#   <out_prefix>_expression.bed.gz and .tbi
#   <out_prefix>_rna_seq_ids.txt
#   <out_prefix>_residualization_summary.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(DESeq2)
  library(edgeR)
  library(Rsamtools)
})

# ------------------------------ CONFIGURATION ---------------------------------

# Input files
counts_path  <- "t2t_aligned_counts_with_sampleIDs.txt"  # gene_id + IID columns
tss_path     <- "hs1.ncbiRefSeq.clean.gene_tss.tsv"
fam_path     <- "T2T_linear/all_chr_T2Tlinear.fam"  # PLINK .fam; IID in column 2
pheno_path   <- "All217SamplePhenotypes_withBatch_t2tPCs_panPCs.txt"
pheno_id_col <- "IID"

# Optional shared keep list (IDs to enforce; NULL to ignore)
shared_keep_path <- NULL

# Optional manual drops by IID
drop_IIDs <- character(0)

# TSS table settings
# Set these column names to match the TSS-table header.
# Supported layouts:
#   A. gene_id + chr + bed_start + bed_end
#   B. gene_id + chr + tss0, where tss0 is already 0-based
#
# Example for option B:
# tss_gene_id_col    <- "gene_id"
# tss_chr_col        <- "chr"
# tss_start_col      <- "tss0"
# tss_end_col        <- NULL
#
# Example for option A:
# tss_gene_id_col    <- "gene_id"
# tss_chr_col        <- "chr"
# tss_start_col      <- "bed_start"
# tss_end_col        <- "bed_end"

tss_gene_id_col <- "gene_id"
tss_chr_col     <- "chr"
tss_start_col   <- "bed_start"
tss_end_col     <- "bed_end"  # Set to NULL when start_col is a 0-based TSS.

# If your TSS file includes transcript-level rows and you want one entry per gene,
# set this to TRUE. It keeps the first TSS encountered per gene after sorting.
deduplicate_tss_by_gene <- FALSE

# Chromosome naming and BED settings
# Genotype chromosomes are expected to use chr1, ..., chr22, chrX, chrY, chrM.
add_chr_prefix_to_tss <- TRUE
standard_chr_levels   <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# Covariate settings
# Sex values are assumed to use PLINK-style coding (1 and 2); value 2 is
# recoded to 1 and all other nonmissing values are recoded to 0.
covar_sex_col <- "sex"
covar_age_col <- "exam_age"

covar_wgs_pc_prefix <- "t2tPC"
covar_wgs_pc_n      <- 10

covar_expr_pc_n <- 10

extra_numeric_covariates <- c()

categorical_covariates <- c(
  "batch"
)

# Analysis and output settings
use_filterByExpr <- TRUE
second_INT <- FALSE

out_prefix      <- "T2T_norm_residualized_counts"
bed_base        <- paste0(out_prefix, "_expression.bed")
ids_out_path    <- paste0(out_prefix, "_rna_seq_ids.txt")
summary_out     <- paste0(out_prefix, "_residualization_summary.txt")

# ------------------------------- FUNCTIONS ------------------------------------

info <- function(...) {
  message("[INFO] ", paste0(..., collapse = ""))
}

warn <- function(...) {
  message("[WARN] ", paste0(..., collapse = ""))
}

abort <- function(...) {
  message("[FAIL] ", paste0(..., collapse = ""))
  stop(call. = FALSE)
}

req_file <- function(x, label = x) {
  if (!file.exists(x)) abort("Missing file: ", label, " at ", x)
}

invnorm <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(r))
}

normalize_chr_names <- function(chr_vec, add_chr_prefix = TRUE) {
  chr_vec <- as.character(chr_vec)

  if (add_chr_prefix) {
    chr_vec <- ifelse(chr_vec == "MT", "chrM", chr_vec)
    chr_vec <- ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
  } else {
    chr_vec <- sub("^chr", "", chr_vec)
    chr_vec <- ifelse(chr_vec == "M", "MT", chr_vec)
  }

  chr_vec
}

build_gene_annot_tss_from_table <- function(
    tss_path,
    gene_id_col,
    chr_col,
    start_col,
    end_col = NULL,
    add_chr_prefix = TRUE,
    keep_chr = NULL,
    dedup_by_gene = TRUE
) {
  dt <- fread(tss_path)

  needed <- c(gene_id_col, chr_col, start_col)
  missing_needed <- setdiff(needed, names(dt))
  if (length(missing_needed)) {
    abort("Missing required TSS columns: ", paste(missing_needed, collapse = ", "))
  }

  if (!is.null(end_col) && !(end_col %in% names(dt))) {
    abort("Requested end_col not found in TSS file: ", end_col)
  }

  out <- dt[, .(
    gene_id_raw = as.character(get(gene_id_col)),
    chr_raw     = as.character(get(chr_col)),
    start_raw   = suppressWarnings(as.integer(get(start_col))),
    end_raw     = if (is.null(end_col)) NA_integer_ else suppressWarnings(as.integer(get(end_col)))
  )]

  # Keep gene IDs exactly as provided in the TSS table.
  out[, gene_id := gene_id_raw]
  out[, chr := normalize_chr_names(chr_raw, add_chr_prefix = add_chr_prefix)]

  # TSS file is already 0-based.
  # If no end column is supplied, assume start_raw is the 0-based TSS coordinate
  # and make a 1-bp BED interval [tss0, tss0+1).
  if (is.null(end_col)) {
    out[, bed_start := start_raw]
    out[, bed_end   := start_raw + 1L]
  } else {
    out[, bed_start := start_raw]
    out[, bed_end   := end_raw]
  }

  out <- out[!is.na(gene_id) & !is.na(chr) & !is.na(bed_start) & !is.na(bed_end)]

  # Basic sanity checks
  out <- out[bed_start >= 0]
  out <- out[bed_end > bed_start]

  if (!is.null(keep_chr)) {
    out <- out[chr %in% keep_chr]
  }

  setorder(out, chr, bed_start, bed_end, gene_id)

  if (dedup_by_gene) {
    dup_n <- sum(duplicated(out$gene_id))
    if (dup_n > 0) {
      warn("TSS table has duplicate gene IDs; keeping first entry per gene (n duplicated = ", dup_n, ").")
      out <- out[!duplicated(gene_id)]
    }
  }

  out <- unique(out[, .(chr, bed_start, bed_end, gene_id)])
  as.data.frame(out)
}

# --------------------------------- WORKFLOW ------------------------------------

req_file(counts_path, "counts")
req_file(tss_path,    "TSS table")
req_file(fam_path,    "FAM")
req_file(pheno_path,  "phenotype")
if (!is.null(shared_keep_path)) req_file(shared_keep_path, "shared_samples.txt")

summary_lines <- c()
add_summary <- function(...) {
  s <- paste0(..., collapse = "")
  summary_lines <<- c(summary_lines, s)
  info(s)
}

# Load count matrix
info("Loading counts ...")
counts_dt <- fread(counts_path)
if (ncol(counts_dt) < 2) abort("Counts need gene_id + ≥1 sample column.")
setnames(counts_dt, 1, "gene_id")
all_iids <- colnames(counts_dt)[-1]
add_summary("Total RNA-seq samples (columns in counts): ", length(all_iids))

# Keep expression gene IDs exactly as provided in the counts file.
counts_dt[, gene_id := as.character(gene_id)]

# Define the analysis sample set in PLINK order
fam <- fread(fam_path, header = FALSE)
if (ncol(fam) < 2) abort("FAM appears malformed; expected ≥2 columns.")
fam_iids <- fam[[2]]
add_summary("Total samples in FAM: ", length(fam_iids))

shared <- NULL
if (!is.null(shared_keep_path)) {
  shared <- read_tsv(shared_keep_path, col_names = FALSE, show_col_types = FALSE)[[1]]
  add_summary("Shared keep list provided (n=", length(shared), ").")
}

keep_iids <- intersect(fam_iids, all_iids)
if (!is.null(shared)) keep_iids <- intersect(keep_iids, shared)
if (length(drop_IIDs)) keep_iids <- setdiff(keep_iids, drop_IIDs)

if (length(keep_iids) == 0) abort("No overlapping samples between counts and FAM after filtering.")
add_summary("WGS-matched samples after drops: ", length(keep_iids))

counts_dt <- counts_dt[, c("gene_id", keep_iids), with = FALSE]

write_tsv(tibble(sample_id = keep_iids), ids_out_path, col_names = FALSE)
info("Wrote sample ID list: ", ids_out_path)

# Filter lowly expressed genes
info("Applying gene filtering (edgeR::filterByExpr) ...")
mat <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(mat) <- counts_dt$gene_id

dge <- DGEList(counts = mat)
if (use_filterByExpr) {
  keep_genes <- filterByExpr(dge)
  add_summary("Genes retained by filterByExpr: ", sum(keep_genes), " / ", nrow(dge))
} else {
  min_counts <- 10
  min_frac   <- 0.2
  nsamp_req  <- ceiling(min_frac * ncol(mat))
  keep_genes <- rowSums(mat >= min_counts) >= nsamp_req
  add_summary("Genes retained by simple count filter: ", sum(keep_genes), " / ", nrow(dge))
}
mat_f <- mat[keep_genes, , drop = FALSE]

# Apply DESeq2 variance-stabilizing transformation
info("Running DESeq2 VST ...")
coldata <- tibble(IID = colnames(mat_f)) %>%
  column_to_rownames("IID")
dds <- DESeqDataSetFromMatrix(countData = round(mat_f), colData = coldata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

add_summary("Genes after filtering (rows in VST matrix): ", nrow(vst_mat))

# Compute expression principal components
expr_pcs_df <- NULL
if (covar_expr_pc_n > 0) {
  info("Computing expression PCs from VST (k=", covar_expr_pc_n, ") ...")
  z <- t(scale(t(vst_mat), center = TRUE, scale = FALSE))
  z[is.na(z)] <- 0
  pc <- prcomp(t(z), center = FALSE, scale. = FALSE)
  k  <- min(covar_expr_pc_n, ncol(pc$x))
  expr_pcs_df <- as_tibble(pc$x[, 1:k, drop = FALSE], rownames = "IID")
  names(expr_pcs_df)[-1] <- paste0("exprPC", 1:k)
  add_summary("Expression PCs computed: ", k)
} else {
  add_summary("Expression PCs disabled (covar_expr_pc_n=0).")
}

# Apply per-gene inverse-normal transformation
info("Applying per-gene inverse-normal transform (INT) ...")
expr_int <- t(apply(vst_mat, 1, invnorm))
rownames(expr_int) <- rownames(vst_mat)
colnames(expr_int) <- colnames(vst_mat)

# Assemble residualization covariates
info("Loading phenotype table for covariates ...")
pheno <- fread(pheno_path) %>% as_tibble()
if (!(pheno_id_col %in% names(pheno))) {
  abort("Phenotype table must have ID column '", pheno_id_col, "'.")
}

cov_df <- tibble(IID = keep_iids)
cov_used_numeric <- c()
cov_used_categorical <- c()

if (!is.null(covar_sex_col)) {
  if (!(covar_sex_col %in% names(pheno))) {
    abort("Sex column '", covar_sex_col, "' not found in phenotype table.")
  }
  tmp <- pheno %>% select(!!sym(pheno_id_col), !!sym(covar_sex_col))
  names(tmp) <- c("IID", "sex_src")
  cov_df <- cov_df %>%
    left_join(tmp, by = "IID") %>%
    mutate(sex = {
      sx <- suppressWarnings(as.integer(as.character(.data$sex_src)))
      ifelse(is.na(sx), NA_integer_, ifelse(sx == 2L, 1L, 0L))
    }) %>%
    select(-sex_src)
  cov_used_numeric <- c(cov_used_numeric, "sex")
} else {
  add_summary("Sex covariate: DISABLED (covar_sex_col=NULL)")
}

if (!is.null(covar_age_col)) {
  if (!(covar_age_col %in% names(pheno))) {
    abort("Age column '", covar_age_col, "' not found in phenotype table.")
  }
  tmp <- pheno %>% select(!!sym(pheno_id_col), !!sym(covar_age_col))
  names(tmp) <- c("IID", "age")
  tmp <- tmp %>% mutate(age = suppressWarnings(as.numeric(age)))
  cov_df <- cov_df %>% left_join(tmp, by = "IID")
  cov_used_numeric <- c(cov_used_numeric, "age")
} else {
  add_summary("Age covariate: DISABLED (covar_age_col=NULL)")
}

if (covar_wgs_pc_n > 0) {
  pc_names <- paste0(covar_wgs_pc_prefix, seq_len(covar_wgs_pc_n))
  missing_pcs <- setdiff(pc_names, names(pheno))
  if (length(missing_pcs)) {
    abort("WGS PC columns missing in phenotype table: ", paste(missing_pcs, collapse = ", "))
  }
  pcs <- pheno %>% select(!!sym(pheno_id_col), all_of(pc_names))
  names(pcs)[1] <- "IID"
  cov_df <- cov_df %>% left_join(pcs, by = "IID")
  cov_used_numeric <- c(cov_used_numeric, pc_names)
} else {
  add_summary("WGS PCs: DISABLED (covar_wgs_pc_n=0)")
}

if (length(extra_numeric_covariates) > 0) {
  missing_extra <- setdiff(extra_numeric_covariates, names(pheno))
  if (length(missing_extra)) {
    abort("Extra numeric covariate columns missing in phenotype table: ",
          paste(missing_extra, collapse = ", "))
  }
  extra <- pheno %>% select(!!sym(pheno_id_col), all_of(extra_numeric_covariates))
  names(extra)[1] <- "IID"
  cov_df <- cov_df %>% left_join(extra, by = "IID")
  cov_used_numeric <- c(cov_used_numeric, extra_numeric_covariates)
}

if (length(categorical_covariates) > 0) {
  missing_cat <- setdiff(categorical_covariates, names(pheno))
  if (length(missing_cat)) {
    abort("Categorical covariate columns missing in phenotype table: ",
          paste(missing_cat, collapse = ", "))
  }
  cat_df <- pheno %>% select(!!sym(pheno_id_col), all_of(categorical_covariates))
  names(cat_df)[1] <- "IID"
  cov_df <- cov_df %>% left_join(cat_df, by = "IID")
  cov_used_categorical <- categorical_covariates
}

if (!is.null(expr_pcs_df)) {
  cov_df <- cov_df %>% left_join(expr_pcs_df, by = "IID")
  cov_used_numeric <- c(cov_used_numeric, names(expr_pcs_df)[-1])
}

add_summary(
  "Numeric covariates requested: ",
  ifelse(length(cov_used_numeric) == 0, "NONE", paste(cov_used_numeric, collapse = ", "))
)
add_summary(
  "Categorical covariates requested: ",
  ifelse(length(cov_used_categorical) == 0, "NONE", paste(cov_used_categorical, collapse = ", "))
)

covariate_df <- cov_df
numeric_cols     <- cov_used_numeric
categorical_cols <- cov_used_categorical

if (length(numeric_cols) > 0) {
  for (nm in numeric_cols) {
    v <- covariate_df[[nm]]
    v <- suppressWarnings(as.numeric(v))
    if (all(is.na(v))) {
      warn("Numeric covariate '", nm, "' is all NA; will be dropped.")
      covariate_df[[nm]] <- NULL
      next
    }
    if (anyNA(v)) {
      m <- mean(v, na.rm = TRUE)
      v[is.na(v)] <- m
    }
    covariate_df[[nm]] <- v
  }

  numeric_cols_present <- intersect(numeric_cols, names(covariate_df))
  if (length(numeric_cols_present) > 0) {
    covariate_df[numeric_cols_present] <- scale(
      covariate_df[numeric_cols_present],
      center = TRUE, scale = FALSE
    )
  }
}

if (length(categorical_cols) > 0) {
  for (nm in categorical_cols) {
    v <- covariate_df[[nm]]
    v_chr <- as.character(v)
    v_chr[is.na(v_chr)] <- "MISSING"
    f <- factor(v_chr)
    if (nlevels(f) <= 1) {
      warn("Categorical covariate '", nm, "' has <=1 level; will be dropped.")
      covariate_df[[nm]] <- NULL
    } else {
      covariate_df[[nm]] <- f
    }
  }
}

rownames(covariate_df) <- covariate_df$IID
covariate_df$IID <- NULL

keep_cols <- sapply(covariate_df, function(v) !all(is.na(v)))
covariate_df <- covariate_df[, keep_cols, drop = FALSE]

if (ncol(covariate_df) == 0) {
  abort("No usable covariates after cleaning; covariate_df is empty.")
}

add_summary(
  "Covariate columns used in model.matrix (after cleaning): ",
  paste(colnames(covariate_df), collapse = ", ")
)

X <- model.matrix(~ ., data = covariate_df)
add_summary("Number of columns in X (including intercept): ", ncol(X))

# Residualize transformed expression
info("Residualizing INT(VST) expression on covariates ...")
Y <- t(expr_int)

XtX <- t(X) %*% X
XtX_inv <- solve(XtX)
P <- X %*% XtX_inv %*% t(X)

Y_fitted <- P %*% Y
Y_resid  <- Y - Y_fitted

expr_resid <- t(Y_resid)
rownames(expr_resid) <- rownames(expr_int)
colnames(expr_resid) <- colnames(expr_int)

if (second_INT) {
  info("Applying second INT to residuals ...")
  expr_resid <- t(apply(expr_resid, 1, invnorm))
  rownames(expr_resid) <- rownames(expr_int)
  colnames(expr_resid) <- colnames(expr_int)
  add_summary("Second INT on residuals: YES")
} else {
  add_summary("Second INT on residuals: NO")
}

# Construct and index the tensorQTL phenotype BED
info("Building gene annotation from TSS table and creating BED ...")
gene_annot <- build_gene_annot_tss_from_table(
  tss_path       = tss_path,
  gene_id_col    = tss_gene_id_col,
  chr_col        = tss_chr_col,
  start_col      = tss_start_col,
  end_col        = tss_end_col,
  add_chr_prefix = add_chr_prefix_to_tss,
  keep_chr       = standard_chr_levels,
  dedup_by_gene  = deduplicate_tss_by_gene
)

if (anyDuplicated(rownames(expr_resid))) {
  dup_ids <- unique(rownames(expr_resid)[duplicated(rownames(expr_resid))])
  abort("Duplicate gene IDs in expr_resid after preprocessing: ",
        paste(head(dup_ids, 20), collapse = ", "))
}

if (anyDuplicated(gene_annot$gene_id)) {
  dup_ids <- unique(gene_annot$gene_id[duplicated(gene_annot$gene_id)])
  abort("Duplicate gene IDs in gene_annot: ",
        paste(head(dup_ids, 20), collapse = ", "))
}

common_genes <- intersect(rownames(expr_resid), gene_annot$gene_id)
if (length(common_genes) < 1000) {
  warn("Few genes after intersect with TSS table; check gene_id formats.")
}
add_summary("Genes in expression & TSS intersection: ", length(common_genes))

expr_use <- expr_resid[common_genes, keep_iids, drop = FALSE]

bed_df <- as_tibble(expr_use, rownames = "gene_id") %>%
  left_join(gene_annot, by = "gene_id") %>%
  filter(!is.na(chr)) %>%
  relocate(chr, bed_start, bed_end, gene_id) %>%
  arrange(factor(chr, levels = standard_chr_levels), bed_start, bed_end)

bed_path <- bed_base
gz_path  <- paste0(bed_base, ".gz")
tbi_path <- paste0(gz_path, ".tbi")

if (file.exists(bed_path)) file.remove(bed_path)
if (file.exists(gz_path))  file.remove(gz_path)
if (file.exists(tbi_path)) file.remove(tbi_path)

header <- paste(c("#chr", "start", "end", "gene_id", keep_iids), collapse = "\t")
writeLines(header, con = bed_path)
fwrite(bed_df, bed_path, sep = "\t", col.names = FALSE, append = TRUE)

Rsamtools::bgzip(bed_path, dest = gz_path, overwrite = TRUE)
Rsamtools::indexTabix(gz_path, format = "bed")

info("Wrote expression BED: ", gz_path, " (+ .tbi)")
add_summary("Final samples in BED: ", length(keep_iids))
add_summary("Final genes in BED: ", nrow(bed_df))

writeLines(summary_lines, con = summary_out)
info("Wrote residualization summary: ", summary_out)
info("DONE.")
