
library(readr)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)
library(DESeq2)
library(ggplot2)
library(GENIE3)
library(tibble)
library(doParallel)
library(doRNG)
library(igraph)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE185263", "file=GSE185263_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000001111111111111111111111111111111111",
               "11111000000000000000000000000000000000000000000000",
               "000000000000000000000000000000000000011111")
sml <- strsplit(gsms, split="")[[1]]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("sepsis","healthy"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:250],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))

########################################################################
filtered_data <- as.matrix(tbl[rownames(tbl) %in% tT$GeneID, ])

sepsis_cols <- rownames(sample_info)[sample_info$Group == "sepsis"]
healthy_cols <- rownames(sample_info)[sample_info$Group == "healthy"]

sepsis_data <- filtered_data[, sepsis_cols, drop = FALSE]
healthy_data <- filtered_data[, healthy_cols, drop = FALSE]
#write.csv(sepsis_data, "/sepsis_data_185263.csv", row.names = TRUE)
#write.csv(healthy_data, "/healthy_data_185263.csv", row.names = TRUE)
#write_csv(tT, "GSE185263_DESeq.csv")
human_tfs <- read_csv('./Data/human_tfs.csv', show_col_types = FALSE)

deseq_info <- tT

get_en_regulators <- function(data, tfs){
  target_genes <- row.names(data)
  is_regulator <- target_genes %in% tfs$Entrez.ID
  regulatoryGenes <- target_genes[is_regulator]
  return(regulatoryGenes)
}

#function that applies GENIE3
get_linkList <- function(expr_data, regulatory_genes){
  weightMat <- GENIE3(expr_data, regulators = regulatory_genes, nCores = 4, verbose = TRUE)
  linklist <- getLinkList(weightMat)
  min_weight <- min(linklist$weight)
  max_weight <- max(linklist$weight)
  linklist$norm_weight <- (linklist$weight - min_weight) / (max_weight - min_weight)
  linklist <- linklist%>%
    dplyr::select(-weight)
  return(linklist)
}

healthy_reg <- get_en_regulators(healthy_data, human_tfs)
sepsis_reg <- get_en_regulators(sepsis_data, human_tfs)

sepsis_linklist <- get_linkList(as.matrix(sepsis_data), sepsis_reg)

sepsis_toplink <- head(sepsis_linklist, 100)

network_graph_sepsis <- graph_from_data_frame(sepsis_toplink, directed = TRUE)
# Open a PNG device
#png("GRN_sepsis_rnaSeq.png", width = 8*300, height = 6*300, res = 300)

# Plot the graph
plot(network_graph_sepsis, vertex.size = 5, 
     vertex.label.cex = 0.8, edge.arrow.size = 0.5, layout = layout_with_fr,
     main = "GRN of key genes in sepsis samples")
healthy_linkList <- get_linkList(as.matrix(healthy_data), healthy_reg)

healthy_toplink <- head(healthy_linkList, 100)
network_graph_healthy <- graph_from_data_frame(healthy_toplink, directed = TRUE)
#png("GRN_healthy_rnaSeq.png", width = 8*300, height = 6*300, res = 300)

# Plot the graph
plot(network_graph_healthy, vertex.size = 5, 
     vertex.label.cex = 0.8, edge.arrow.size = 0.5, layout = layout_with_fr,
     main = "GRN of key genes for sepsis in healthy samples")

diff_info <- tT
diff_data <- diff_info%>%
  dplyr::select(c(Symbol, log2FoldChange, padj))%>%
  dplyr::filter(padj <= 0.05 & Symbol != "NA")%>%
  dplyr::select(-padj)

diff_data <- diff_data %>%
  column_to_rownames(var = "Symbol")

geneList <- diff_data$log2FoldChange
names(geneList) <- row.names(diff_data)
gene_symbol <- rownames(diff_data)


entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbol,
                     column = "ENTREZID", keytype = "SYMBOL",
                     multiVals = "first")

result <- enrichGO(gene = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

barplot(result, title = "GO enrichment of Molecular Function")

#change to other biological terms for other GO enrichment analysis like "BP", "CC"
