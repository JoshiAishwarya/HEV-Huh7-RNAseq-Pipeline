library(tximport)
library(DESeq2)
library(pheatmap)
library(matrixStats)
library(org.Hs.eg.db)
library(apeglm)

# 1. Setup Output Directory Boundaries on Windows D Drive
out_dir <- "D:/HEV_Analysis/Human/Results/DESeq2_Results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# 2. Reconstruct the Metadata Sheet exactly as shown in your image
samples <- data.frame(
  sample = c("SRR5407646", "SRR5407647", "SRR5407648", "SRR5407653", "SRR5407654", "SRR5407655", "SRR5407656"),
  time = as.factor(c("72h", "72h", "72h", "24h", "24h", "24h", "24h")),
  condition = as.factor(c("infected", "infected", "control", "control", "control", "infected", "infected")),
  stringsAsFactors = FALSE
)

# Synchronize row order alphabetically to match your hard drive files exactly
samples <- samples[order(samples$sample), ]

# 3. Reference Path Pointers Straight to Nextflow Outputs
files <- file.path("D:/HEV_Analysis/Human/Results/salmon_quant", samples$sample, "quant.sf")
names(files) <- samples$sample

# 4. Extract Identifiers and Build Transcript Mapper (Single File Vector Index Fix)
first_quant <- read.delim(files[1], stringsAsFactors=FALSE)
clean_tx_names <- gsub("\\..*", "", first_quant$Name)
tx2gene <- data.frame(TXNAME=clean_tx_names, GENEID=clean_tx_names)

# 5. Ingest Clean Expression Abundance Values
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)

# 6. Set Reference Levels (Control and 24h are the baselines)
samples$condition <- relevel(samples$condition, ref="control")
samples$time <- relevel(samples$time, ref="24h")

# ==============================================================================
# MASTER MODEL: INTERACTION DESIGN (Uses all 7 samples to fix zero-variance)
# ==============================================================================
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~ time + condition + time:condition)

# Filter out low-count background noise genes upfront to stabilize shrinkage
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

# Get the available model names to ensure correct matching coefficients
res_names <- resultsNames(dds)

# ------------------------------------------------------------------------------
# EXTRACTION 1: PURE 24 HOUR PHENOTYPE DIFFERENCES
# ------------------------------------------------------------------------------
res_24h_shrunk <- lfcShrink(dds, coef="condition_infected_vs_control", type="apeglm")
df_24h <- as.data.frame(res_24h_shrunk)
df_24h$transcript_id <- gsub("\\..*", "", rownames(df_24h))
df_24h$gene_symbol <- mapIds(org.Hs.eg.db, keys=df_24h$transcript_id, column="SYMBOL", keytype="REFSEQ", multiVals="first")
df_24h <- df_24h[, c("gene_symbol", "transcript_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write.csv(df_24h, file.path(out_dir, "overall_master_results.csv"), row.names=FALSE)
sig_24h <- df_24h[!is.na(df_24h$gene_symbol) & !is.na(df_24h$padj) & df_24h$padj < 0.05, ]
write.csv(head(sig_24h[order(sig_24h$log2FoldChange, decreasing = TRUE), ], 20), file.path(out_dir, "time24h_top20_upregulated.csv"), row.names=FALSE)
write.csv(head(sig_24h[order(sig_24h$log2FoldChange, decreasing = FALSE), ], 20), file.path(out_dir, "time24h_top20_downregulated.csv"), row.names=FALSE)

# ------------------------------------------------------------------------------
# EXTRACTION 2: PURE 72 HOUR PHENOTYPE DIFFERENCES (Calculated via interaction coefficient)
# ------------------------------------------------------------------------------
res_72h_shrunk <- lfcShrink(dds, coef="time72h.conditioninfected", type="apeglm")
df_72h <- as.data.frame(res_72h_shrunk)
df_72h$transcript_id <- gsub("\\..*", "", rownames(df_72h))
df_72h$gene_symbol <- mapIds(org.Hs.eg.db, keys=df_72h$transcript_id, column="SYMBOL", keytype="REFSEQ", multiVals="first")
df_72h <- df_72h[, c("gene_symbol", "transcript_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write.csv(df_72h, file.path(out_dir, "deseq2_72h_results.csv"), row.names=FALSE)
sig_72h <- df_72h[!is.na(df_72h$gene_symbol) & !is.na(df_72h$padj) & df_72h$padj < 0.05, ]
write.csv(head(sig_72h[order(sig_72h$log2FoldChange, decreasing = TRUE), ], 20), file.path(out_dir, "time72h_top20_upregulated.csv"), row.names=FALSE)
write.csv(head(sig_72h[order(sig_72h$log2FoldChange, decreasing = FALSE), ], 20), file.path(out_dir, "time72h_top20_downregulated.csv"), row.names=FALSE)

# ------------------------------------------------------------------------------
# EXTRACTION 3: THE COMPREHENSIVE OVERALL RESULTS COMBINATION
# ------------------------------------------------------------------------------
res_overall <- results(dds, name="condition_infected_vs_control")
df_overall <- as.data.frame(res_overall)
df_overall$transcript_id <- gsub("\\..*", "", rownames(df_overall))
df_overall$gene_symbol <- mapIds(org.Hs.eg.db, keys=df_overall$transcript_id, column="SYMBOL", keytype="REFSEQ", multiVals="first")
df_overall <- df_overall[, c("gene_symbol", "transcript_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
write.csv(df_overall, file.path(out_dir, "deseq2_overall_results.csv"), row.names=FALSE)
sig_overall <- df_overall[!is.na(df_overall$gene_symbol) & !is.na(df_overall$padj) & df_overall$padj < 0.05, ]
write.csv(head(sig_overall[order(sig_overall$log2FoldChange, decreasing = TRUE), ], 20), file.path(out_dir, "overall_top20_upregulated.csv"), row.names=FALSE)
write.csv(head(sig_overall[order(sig_overall$log2FoldChange, decreasing = FALSE), ], 20), file.path(out_dir, "overall_top20_downregulated.csv"), row.names=FALSE)

# ==============================================================================
# GRAPHICS ENGINE GENERATION
# ==============================================================================
vsd <- vst(dds, blind=FALSE)

png(file.path(out_dir, "pca_plot.png"), width=900, height=800, res=150)
print(plotPCA(vsd, intgroup=c("condition", "time")))
dev.off()

png(file.path(out_dir, "ma_plot.png"), width=900, height=800, res=150)
plotMA(res_overall, ylim=c(-5,5), main="MA Plot: Infected vs Control")
dev.off()

top_genes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing=TRUE), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
plot_gene_names <- mapIds(org.Hs.eg.db, keys=rownames(mat), column="SYMBOL", keytype="REFSEQ", multiVals="first")
plot_gene_names[is.na(plot_gene_names)] <- rownames(mat)[is.na(plot_gene_names)] 
rownames(mat) <- plot_gene_names
df_anno <- as.data.frame(colData(dds)[,c("condition","time")])
png(file.path(out_dir, "expression_heatmap.png"), width=1000, height=1200, res=150)
pheatmap(mat, annotation_col=df_anno, show_rownames=TRUE, main="Top 50 Most Variable Genes")
dev.off()

png(file.path(out_dir, "volcano_24h.png"), width=800, height=800, res=150)
plot(res_24h_shrunk$log2FoldChange, -log10(res_24h_shrunk$pvalue), pch=20, main="DE Analysis: HEV at 24 hrs", 
     xlab="log2 Fold Change", ylab="-log10 p-value", col=ifelse(abs(res_24h_shrunk$log2FoldChange) > 1 & res_24h_shrunk$padj < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1, 1), col="blue", lty=2)
dev.off()

png(file.path(out_dir, "volcano_72h.png"), width=800, height=800, res=150)
plot(res_72h_shrunk$log2FoldChange, -log10(res_72h_shrunk$pvalue), pch=20, main="DE Analysis: HEV at 72 hrs", 
     xlab="log2 Fold Change", ylab="-log10 p-value", col=ifelse(abs(res_72h_shrunk$log2FoldChange) > 1 & res_72h_shrunk$padj < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1, 1), col="blue", lty=2)
dev.off()

cat("\nSUCCESS: All files generated using the unified interaction baseline model!\n")


#==============================================================================
# 2. Safely read and rename the 24h master file to its correct uniform title
if (file.exists(file.path(out_dir, "overall_master_results.csv"))) {
  file.rename(from = file.path(out_dir, "overall_master_results.csv"),
              to = file.path(out_dir, "deseq2_24h_results.csv"))
  message("SUCCESS: File renamed cleanly to deseq2_24h_results.csv!")
} else {
  message("Notice: File was already adjusted or could not be located.")
}
#==============================================================================