#!/usr/bin/env Rscript

# Real-ish DESeq2 analysis for the rare_disease_multiomics_reporter project,
# using gene-wise dispersion estimates only (simpler for tiny toy datasets).
#
# Usage:
#   Rscript 01_deseq2_analysis.R <counts.tsv> <output_dir>
#
# Assumptions:
#   - counts.tsv: genes in rows, samples in columns
#   - at least four columns of counts
#   - first two columns = controls, next two columns = treated

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 01_deseq2_analysis.R <counts.tsv> <output_dir>", call. = FALSE)
}

counts_path <- args[[1]]
out_dir     <- args[[2]]

cat("=== DESeq2 analysis (toy, gene-wise dispersions) ===\n")
cat("Counts file:", counts_path, "\n")
cat("Output dir :", out_dir, "\n")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Read counts
counts <- read.table(
  counts_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

if (ncol(counts) < 4) {
  stop("Expected at least four columns of counts (2 control, 2 treated).", call. = FALSE)
}

sample_names <- colnames(counts)

# First two = control, next two = treated
condition <- factor(c("control", "control", "treated", "treated"),
                    levels = c("control", "treated"))
coldata <- data.frame(
  row.names = sample_names[1:4],
  condition = condition
)

dds <- DESeqDataSetFromMatrix(
  countData = counts[, 1:4],
  colData   = coldata,
  design    = ~ condition
)

# Use the simpler, gene-wise-only workflow recommended in the error message
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

res <- results(dds, contrast = c("condition", "treated", "control"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, "deseq2_results.tsv")
write.table(res_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote DESeq2 results to:", out_file, "\n")

# Volcano plot
volcano_file <- file.path(out_dir, "volcano_deseq2.png")
png(volcano_file, width = 800, height = 600)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.7) +
  xlab("log2 fold change (treated vs control)") +
  ylab("-log10(p-value)") +
  ggtitle("DESeq2 volcano plot (toy, gene-wise dispersions)") +
  theme_minimal()
dev.off()
cat("Saved volcano plot to:", volcano_file, "\n")

cat("DESeq2 analysis complete.\n")

