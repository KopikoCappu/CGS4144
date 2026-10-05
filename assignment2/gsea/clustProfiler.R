library(clusterProfiler)
library(org.At.tair.db)

df <- read.csv("~/CGS4144/Project/CGS4144/results/gsea_ranked_genes.csv")

gene_list <- df$Effect
names(gene_list) <- df$Gene

gene_list <- sort(gene_list, decreasing = TRUE)

gsea <- gseGO(
  geneList = gene_list,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pvalueCutoff = 0.05
)

results <- as.data.frame(gsea)

write.csv(
  results,
  "clusterProfiler_results.csv",
  row.names = FALSE
)

head(results)