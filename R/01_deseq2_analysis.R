#!/usr/bin/env Rscript

# Simple placeholder expression "analysis" for the rare_disease_multiomics_reporter project.
# This is NOT real DESeq2, just a light-weight demonstration of:
#  - reading a counts matrix
#  - doing a basic transformation
#  - writing results to results/expression/

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 01_deseq2_analysis.R <counts.tsv> <output_dir>", call. = FALSE)
}

counts_path <- args[[1]]
out_dir     <- args[[2]]

cat("=== R expression analysis placeholder ===\n")
cat("Counts file:", counts_path, "\n")
cat("Output dir :", out_dir, "\n")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Read counts: genes in rows, samples in columns
counts <- read.table(counts_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

if (ncol(counts) < 2) {
  stop("Expected at least two columns of counts (e.g. sample1, sample2).", call. = FALSE)
}

# Very simple "log2 fold change" between column 2 and column 1
log2fc <- log2((counts[, 2] + 1) / (counts[, 1] + 1))

res <- data.frame(
  gene    = rownames(counts),
  baseMean = rowMeans(counts),
  log2FC   = as.numeric(log2fc),
  pvalue   = 1.0,
  padj     = 1.0,
  row.names = NULL
)

out_file <- file.path(out_dir, "expression_results.tsv")
write.table(res, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote expression results to:", out_file, "\n")

