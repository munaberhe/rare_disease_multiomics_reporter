#!/usr/bin/env Rscript

# GO/Pathway enrichment using clusterProfiler for the rare_disease_multiomics_reporter project.
# Usage:
#   Rscript 02_pathway_analysis.R <deseq2_results.tsv> <output_dir>

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 02_pathway_analysis.R <deseq2_results.tsv> <output_dir>", call. = FALSE)
}

res_path <- args[[1]]
out_dir  <- args[[2]]

cat("=== GO enrichment (toy) ===\n")
cat("DESeq2 results file:", res_path, "\n")
cat("Output dir          :", out_dir, "\n")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

res <- read.table(res_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter for "significant" DE genes (toy thresholds)
sig <- subset(res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

cat("Number of significant genes (padj < 0.05 & |log2FC| > 1):", nrow(sig), "\n")

if (nrow(sig) == 0) {
  cat("No significant genes found; skipping GO enrichment.\n")
  quit(save = "no")
}

genes <- unique(sig$gene)

# Try mapping from SYMBOL or ENSEMBL-like IDs to ENTREZID
mapping <- tryCatch(
  bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  ),
  error = function(e) {
    cat("Error in bitr():", conditionMessage(e), "\n")
    NULL
  }
)

if (is.null(mapping) || nrow(mapping) == 0) {
  cat("No genes could be mapped to ENTREZ IDs; skipping GO enrichment.\n")
  quit(save = "no")
}

cat("Mapped", nrow(mapping), "genes to ENTREZ IDs.\n")

ego <- enrichGO(
  gene          = unique(mapping$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  cat("GO enrichment returned no significant terms.\n")
  quit(save = "no")
}

ego_df <- as.data.frame(ego)
go_file <- file.path(out_dir, "go_enrichment_BP.tsv")
write.table(ego_df, file = go_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote GO BP enrichment results to:", go_file, "\n")

# Simple barplot of top terms
barplot_file <- file.path(out_dir, "go_bp_barplot.png")
png(barplot_file, width = 900, height = 600)
barplot(
  ego,
  showCategory = 10,
  title = "Top GO BP terms (toy)",
  font.size = 10
)
dev.off()
cat("Saved GO BP barplot to:", barplot_file, "\n")

cat("GO enrichment complete.\n")

