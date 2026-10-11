library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler) 

# 1. Establish your directory destination lines
out_dir <- "D:/HEV_Analysis/Human/Results/DESeq2_Results"

# ==============================================================================
# 2. GLOBAL ENRICHMENT PROFILE RUN: 24-HOUR ACUTE WINDOW
# ==============================================================================
master_24h <- read.csv(file.path(out_dir, "deseq2_24h_results.csv"))

sig_up_24h <- master_24h[!is.na(master_24h$gene_symbol) & 
                           !is.na(master_24h$padj) & 
                           master_24h$padj < 0.05 & 
                           master_24h$log2FoldChange > 0, ]

master_genes_24h <- as.character(sig_up_24h$gene_symbol)

# FIXED: Removed the ontology parameter completely to leverage the engine default
go_master_24h <- clusterProfiler::enrichGO(gene          = master_genes_24h,
                                           OrgDb         = org.Hs.eg.db,
                                           keyType       = 'SYMBOL',
                                           pAdjustMethod = 'BH',
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05)

# Save the comprehensive 24h Dotplot
png(file.path(out_dir, "master_enrichment_dotplot_24h.png"), width=1400, height=1000, res=150)
print(dotplot(go_master_24h, showCategory=15, title="Global Host Pathways Activated at 24h (Full Dataset)"))
dev.off()


# ==============================================================================
# 3. GLOBAL ENRICHMENT PROFILE RUN: 72-HOUR TERMINAL PHASE
# ==============================================================================
master_72h <- read.csv(file.path(out_dir, "deseq2_72h_results.csv"))

sig_up_72h <- master_72h[!is.na(master_72h$gene_symbol) & 
                           !is.na(master_72h$padj) & 
                           master_72h$padj < 0.05 & 
                           master_72h$log2FoldChange > 0, ]

master_genes_72h <- as.character(sig_up_72h$gene_symbol)

# FIXED: Removed the ontology parameter completely to leverage the engine default
go_master_72h <- clusterProfiler::enrichGO(gene          = master_genes_72h,
                                           OrgDb         = org.Hs.eg.db,
                                           keyType       = 'SYMBOL',
                                           pAdjustMethod = 'BH',
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05)

# Save the comprehensive 72h Dotplot
png(file.path(out_dir, "master_enrichment_dotplot_72h.png"), width=1400, height=1000, res=150)
print(dotplot(go_master_72h, showCategory=15, title="Global Host Pathways Activated at 72h (Full Dataset)"))
dev.off()

message("SUCCESS: Both master pathway charts compiled cleanly on your D drive!")
