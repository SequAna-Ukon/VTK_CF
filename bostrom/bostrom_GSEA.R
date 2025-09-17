library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db) 


sig_genes_unique <- read.table("sig_genes_unique.txt", stringsAsFactors = FALSE)[,1]
sig_genes_redundant <- read.table("sig_genes_redundant.txt", stringsAsFactors = FALSE)[,1]


# GO enrichment

gene.df.unique <- bitr(sig_genes_unique,
                       fromType = "ENSEMBL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
ego_unique <- enrichGO(gene         = gene.df.unique$ENTREZID,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",
                       pAdjustMethod= "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ego_unique), "GO_enrichment_unique.csv", row.names = FALSE)

# Plot
pdf("GO_barplot_unique.pdf")
barplot(ego_unique, showCategory = 20, title = "GO Enrichment (Unique)")
dev.off()

pdf("GO_dotplot_unique.pdf")
dotplot(ego_unique, showCategory = 20, title = "GO Enrichment (Unique)")
dev.off()

gene.df.redundant <- bitr(sig_genes_redundant,
                          fromType = "ENSEMBL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)
ego_redundant <- enrichGO(gene         = gene.df.redundant$ENTREZID,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "BP",
                          pAdjustMethod= "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ego_redundant), "GO_enrichment_redundant.csv", row.names = FALSE)

# Plot
pdf("GO_barplot_redundant.pdf")
barplot(ego_redundant, showCategory = 20, title = "GO Enrichment (Redundant)")
dev.off()

pdf("GO_dotplot_redundant.pdf")
dotplot(ego_redundant, showCategory = 20, title = "GO Enrichment (Redundant)")
dev.off()



# Run KEGG enrichment for unique genes
ekegg_unique <- enrichKEGG(gene         = gene.df.unique$ENTREZID,
                           organism     = 'hsa',     # human
                           keyType      = "kegg",
                           pvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ekegg_unique), "KEGG_enrichment_unique.csv", row.names = FALSE)

# Plot
pdf("KEGG_barplot_unique.pdf")
barplot(ekegg_unique, showCategory = 20, title = "KEGG Enrichment (Unique)")
dev.off()

pdf("KEGG_dotplot_unique.pdf")
dotplot(ekegg_unique, showCategory = 20, title = "KEGG Enrichment (Unique)")
dev.off()

# Run KEGG enrichment for redundant genes
ekegg_redundant <- enrichKEGG(gene         = gene.df.redundant$ENTREZID,
                              organism     = 'hsa',
                              keyType      = "kegg",
                              pvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ekegg_redundant), "KEGG_enrichment_redundant.csv", row.names = FALSE)

# Plot
pdf("KEGG_barplot_redundant.pdf")
barplot(ekegg_redundant, showCategory = 20, title = "KEGG Enrichment (Redundant)")
dev.off()

pdf("KEGG_dotplot_redundant.pdf")
dotplot(ekegg_redundant, showCategory = 20, title = "KEGG Enrichment (Redundant)")
dev.off()
