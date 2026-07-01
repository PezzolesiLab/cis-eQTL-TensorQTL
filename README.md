# cis-eQTL-TensorQTL

# 01_prepare_tensorqtl_expression_bed_chm13.R - Prepare residualized RNA-seq expression phenotypes for tensorQTL (CHM13/hs1)
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

# 02_run_tensorqtl_mapping.py - Run cis-eQTL mapping with tensorQTL
 # Analyses performed
# 1. Gene-level cis-QTL permutation mapping.
# 2. Nominal cis-QTL mapping, with tensorQTL parquet output by chromosome.
# 3. Conditional cis-QTL mapping to identify independent signals.
# 4. Per-chromosome and overall run summaries.
