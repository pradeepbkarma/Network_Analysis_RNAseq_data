library(tidyverse)
library(readr)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

des_seq_file <- read_csv("./Des_seq_GSE5772.csv", show_col_types = FALSE)
sepsis_diffExpr <- des_seq_file%>%
  filter(adj.P.Val <= 0.05)%>%
  dplyr::select(c(Gene.symbol, logFC))

sepsis_diffExpr <- sepsis_diffExpr%>%
  filter(Gene.symbol != "NA")%>%
  column_to_rownames(var = "Gene.symbol")

geneList <- sepsis_diffExpr$logFC
names(geneList) <- row.names(sepsis_diffExpr)
gene_symbol <- rownames(sepsis_diffExpr)


entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbol,
                     column = "ENTREZID", keytype = "SYMBOL",
                     multiVals = "first")
valid_entrez_id <- entrez_ids[!is.na(entrez_ids)]
geneList <- geneList[names(geneList) %in% names(valid_entrez_id)]
names(geneList) <- valid_entrez_id

result <- enrichGO(gene = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)


barplot(result, title = "GO enrichment of Biological Processes")

barplot(result, title = "GO enrichment of Cellular Components")


