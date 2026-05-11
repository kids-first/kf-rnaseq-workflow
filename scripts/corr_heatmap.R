#!/usr/bin/env Rscript
#'
#' Sample–sample correlation analysis using RSEM TPM expression
#'
#' This script computes and visualizes sample–sample correlations from
#' gene-level TPM expression values produced by RSEM. Expression values
#' are log2-transformed, filtered to retain expressed genes, and used to
#' calculate Pearson correlations between samples.
#'
#' The resulting correlation matrix is clustered using correlation
#' distance (1 − Pearson correlation) and visualized as a heatmap.
#' The script is intended for exploratory quality control and comparison
#' of two sets of RNA-seq samples.
#'
#' Main steps:
#'   1. Read RSEM gene-level TPM files (TSV or TSV.GZ).
#'   2. Normalize gene identifiers by removing version suffixes.
#'   3. Subset to a provided gene list.
#'   4. Filter low-expression genes.
#'   5. Log2-transform TPM values.
#'   6. Compute sample–sample Pearson correlations.
#'   7. Cluster samples using correlation distance and generate a heatmap.
#'
#' Notes:
#'   - Correlation-based similarity is robust to global expression shifts
#'     and is appropriate for cross-sample comparisons.
#'   - Clustering uses 1 − correlation as a distance metric to ensure
#'     consistency between similarity and hierarchical clustering.
#'
#' Usage:
#'   sample_correlation_heatmap.R \\
#'     --gene_list_path genes.tsv \\
#'     --setA_name SetA --files_setA A1.tsv.gz A2.tsv.gz \\
#'     --setB_name SetB --files_setB B1.tsv.gz B2.tsv.gz \\
#'     --output_basename SetA_vs_SetB
#'
#' Outputs:
#'   - <output_basename>.sample_correlation_heatmap.csv
#'   - <output_basename>.sample_correlation_heatmap.pdf
#'   - <output_basename>.removed_genes_tpm.csv
#'
#' Requirements:
#'   R >= 4.0
#'   Packages: readr, dplyr, tidyr, purrr, pheatmap, argparse
#'

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(pheatmap)
  library(argparse)
})

#' Read RSEM TPM file (optionally gzipped)
#'
#' @param path Path to a TSV or TSV.GZ file with columns `gene_id` and `TPM`
#' @return A tibble with columns `gene_id` and `TPM`
read_rsem_tpm_gz <- function(path) {
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    select(gene_id, TPM)
}

# ---- Args ----
p <- ArgumentParser()
p$add_argument("--gene_list_path", required = TRUE)

p$add_argument("--setA_name", default = "setA")
p$add_argument("--files_setA", nargs = "+", required = TRUE)

p$add_argument("--setB_name", default = "setB")
p$add_argument("--files_setB", nargs = "+", required = TRUE)

p$add_argument("--output_basename", default = NULL,
               help = "Basename for output files (no extension).")

args <- p$parse_args()

gene_list_path <- args$gene_list_path
setA_name <- args$setA_name
setB_name <- args$setB_name
out_basename <- args$output_basename
if (is.null(out_basename) || !nzchar(out_basename)) {
  out_basename <- paste0(setA_name, "v", setB_name)
}

files_setA <- setNames(args$files_setA, paste0("A", seq_along(args$files_setA)))
files_setB <- setNames(args$files_setB, paste0("B", seq_along(args$files_setB)))

files <- c(files_setA, files_setB)
# ---- END Args ----


# Long -> wide using a versionless gene key
tpm_long <- imap_dfr(files, ~ read_rsem_tpm_gz(.x) %>% mutate(sample = .y)) %>%
  mutate(
    gene_ensg = sub("_.*$", "", gene_id),        # ENSG...<version>
    gene_ensg = sub("\\..*$", "", gene_ensg),    # ENSG... (no version)
  ) %>%
  select(gene_ensg, sample, TPM)

# If the same gene_ensg appears more than once per sample, collapse it (usually not needed, but safe)
tpm_long <- tpm_long %>%
  group_by(gene_ensg, sample) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

tpm_mat <- tpm_long %>%
  pivot_wider(names_from = sample, values_from = TPM, values_fill = 0) %>%
  column_to_rownames("gene_ensg") %>%
  as.matrix()

genes_tbl <- readr::read_tsv(
  gene_list_path,
  col_names = c("gene_id", "symbol"),
  show_col_types = FALSE
)

# Keep list as Ensembl-without-version (robust)
genes_keep <- unique(sub("\\..*$", "", genes_tbl$gene_id))

# Extract Ensembl-without-version from tpm_mat rownames like "ENSG....15_TSPAN6"
tpm_ensg <- sub("_.*$", "", rownames(tpm_mat))   # "ENSG....15"
tpm_ensg <- sub("\\..*$", "", tpm_ensg)          # "ENSG...."

# Subset
tpm_mat <- tpm_mat[tpm_ensg %in% genes_keep, , drop = FALSE]

# Optional: filter low-expression genes (recommended for correlations)
keep <- rowSums(tpm_mat >= 1) >= 2

removed_rows <- tpm_mat[!keep, , drop = FALSE]
write.csv(
  removed_rows,
  paste0(out_basename, ".removed_genes_tpm.csv"),
  quote = FALSE
)

tpm_mat <- tpm_mat[keep, , drop = FALSE]

# Samples x genes, log-transform
x <- t(log2(tpm_mat + 1))
dim(x)

# Sample–sample correlation
system.time(cors <- cor(t(x), method = "pearson", use = "pairwise.complete.obs"))
write.csv(
  cors,
  paste0(out_basename, ".sample_correlation_heatmap.csv"),
  quote = FALSE
)

# Optional annotation (edit the rule to match your sample names)
ann <- data.frame(Set = ifelse(grepl("^A", rownames(cors)), setA_name, setB_name))
rownames(ann) <- rownames(cors)

# ---- CORRELATION DISTANCE FOR CLUSTERING ----
dist_corr <- as.dist(1 - cors)

# Hierarchical clustering can vary slightly due to ties. Set seed for reproducibility
set.seed(1)
# Pearson correlation on log2(TPM+1); clustering uses 1 - correlation distance
pheatmap(
  cors,
  clustering_distance_rows = dist_corr,
  clustering_distance_cols = dist_corr,
  annotation_col = ann,
  annotation_row = ann,
  main = "Sample correlation (Pearson) on log2(TPM+1)",
  filename = paste0(out_basename, ".sample_correlation_heatmap.pdf"),
  width = 8,
  height = 8
)

