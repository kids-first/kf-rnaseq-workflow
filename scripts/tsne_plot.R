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

# final safety: remove any remaining non-finite values
ok <- is.finite(x)
if (any(!ok)) {
  x <- x[, colSums(!ok) == 0, drop = FALSE]
}

sum(is.na(tpm_mat))
sum(!is.finite(tpm_mat))
dim(tpm_mat)

# t-SNE (perplexity must be small for ~10 samples)
set.seed(1)
perp <- max(2, min(30, floor((nrow(x) - 1) / 3)))

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
