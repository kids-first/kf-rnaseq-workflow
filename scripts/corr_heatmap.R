suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(pheatmap)
  library(argparse)
})

read_rsem_tpm_gz <- function(path) {
  # Works for .gz directly
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
write.csv(removed_rows, "removed_genes_tpm.csv", quote = FALSE)

tpm_mat <- tpm_mat[keep, , drop = FALSE]

# Samples x genes, log-transform
x <- t(log2(tpm_mat + 1))
dim(x)

# Sample–sample correlation
system.time(cors <- cor(t(x), method = "pearson", use = "pairwise.complete.obs"))
write.csv(cors, paste0(out_basename, ".sample_correlation_heatmap.csv"), quote = FALSE)
dim(cors)

# Optional annotation (edit the rule to match your sample names)
ann <- data.frame(Set = ifelse(grepl("^A", rownames(cors)), setA_name, setB_name))
rownames(ann) <- rownames(cors)

pheatmap(
  cors,
  annotation_col = ann,
  annotation_row = ann,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Sample correlation (Pearson) on log2(TPM+1)",
  filename = paste0(out_basename, ".sample_correlation_heatmap.pdf"),
  width = 8,
  height = 8
)


