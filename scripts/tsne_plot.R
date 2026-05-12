#!/usr/bin/env Rscript
#'
#' t-SNE visualization of RNA-seq samples using RSEM TPM values
#'
#' This script performs an exploratory t-SNE analysis of RNA-seq samples
#' based on gene-level TPM expression values produced by RSEM. Expression
#' values are log2-transformed, filtered to retain expressed genes, scaled
#' per gene, and embedded using t-SNE with Euclidean distance.
#'
#' The script is designed for sample-level quality control and comparison
#' between two sets of samples. It supports gzip-compressed input files,
#' robust gene ID normalization, and produces a publication-ready t-SNE
#' scatter plot.
#'
#' Main steps:
#'   1. Read RSEM gene-level TPM files (TSV or TSV.GZ).
#'   2. Normalize gene identifiers by removing version suffixes.
#'   3. Filter low-expression and zero-variance genes.
#'   4. Log2-transform and scale expression values per gene.
#'   5. Compute a t-SNE embedding of samples.
#'   6. Generate a labeled t-SNE plot colored by sample set.
#'
#' Notes:
#'   - t-SNE is an exploratory visualization; distances and clusters should
#'     not be over-interpreted, especially for small sample sizes.
#'   - Euclidean distance is used on scaled log2(TPM+1) values.
#'
#' Usage:
#'   tsne_samples.R \\
#'     --gene_list_path genes.tsv \\
#'     --setA_name SetA --files_setA A1.tsv.gz A2.tsv.gz \\
#'     --setB_name SetB --files_setB B1.tsv.gz B2.tsv.gz \\
#'     --output_basename SetA_vs_SetB
#'
#' Outputs:
#'   - <output_basename>.tsne_plot.png
#'
#' Requirements:
#'   R >= 4.0, packages: readr, dplyr, tidyr, purrr, Rtsne, ggplot2, ggrepel
#'

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(Rtsne)
  library(ggplot2)
  library(ggrepel)
  library(argparse)
})


#' Read an RSEM TPM file (optionally gzipped)
#'
#' Reads an RSEM gene-level expression file containing TPM values.
#' The file may be plain TSV or gzip-compressed (.gz). Only the
#' `gene_id` and `TPM` columns are retained.
#'
#' @param path Character string giving the path to a TSV or TSV.GZ file
#'   with columns `gene_id` and `TPM`.
#'
#' @return A tibble with two columns:
#' \describe{
#'   \item{gene_id}{Gene identifier (e.g., Ensembl ID)}
#'   \item{TPM}{Transcripts Per Million values}
#' }
#'
#' @examples
#' \dontrun{
#' tpm <- read_rsem_tpm_gz("sample.rsem.genes.results.gz")
#' head(tpm)
#' }
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

# Filter low-expression genes (tune if desired)
keep <- rowSums(tpm_mat >= 1) >= 2
tpm_mat <- tpm_mat[keep, , drop = FALSE]

# Samples x genes, log-transform
x <- t(log2(tpm_mat + 1))

# remove zero-variance genes BEFORE scaling (prevents NaN)
nzv <- apply(x, 2, sd) > 0
x <- x[, nzv, drop = FALSE]

# scale genes
x <- scale(x)

# Warn for small groups
if (nrow(x) < 6) {
  warning("Very small sample size for t-SNE; interpretation may be unstable.")
}


# final safety: remove any remaining non-finite values
ok <- is.finite(x)
if (any(!ok)) {
  x <- x[, colSums(!ok) == 0, drop = FALSE]
}

# t-SNE (perplexity must be small for ~10 samples)
set.seed(1)
perp <- max(2, min(30, floor((nrow(x) - 1) / 3)))

# t-SNE uses Euclidean distance on scaled log2(TPM+1)
if (any(duplicated(x))) {
  warning("Duplicate samples detected in expression matrix.")
}
ts <- Rtsne(x, perplexity = perp, pca = TRUE, check_duplicates = FALSE)

emb <- as.data.frame(ts$Y)
colnames(emb) <- c("tSNE1", "tSNE2")
emb$sample <- rownames(x)
emb$set <- ifelse(emb$sample %in% names(files_setA), setA_name, setB_name)


plt <- ggplot(emb, aes(tSNE1, tSNE2, color = set, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  theme_classic() +
  labs(title = "t-SNE of samples using RSEM TPM",
       subtitle = "log2(TPM+1), scaled per gene")

ggsave(paste0(out_basename, ".tsne_plot.png"),
       plot = plt, width = 7, height = 5, dpi = 300)
