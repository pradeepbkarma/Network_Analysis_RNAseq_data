
library(readr)
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)

#### load data and reformat data ####
data <- read_csv('./Data/GDS2808.soft.csv', show_col_types = FALSE)
exp_data <- data[!duplicated(data$IDENTIFIER),]
exp_data <- exp_data%>%
  dplyr::select(-'ID_REF')%>%
  column_to_rownames(var = 'IDENTIFIER')

#### data preprocessing and filtering ####
handle_na_values <- function (data){
  na_counts_per_row <- rowSums(is.na(data))
  threshold<- ncol(data)/2
  
  unwanted_rows <- sum(na_counts_per_row > threshold)
  print(paste("There are ", unwanted_rows, "rows with more than 50% missing values."))
  exp_data_clean <- data %>%
    filter(rowSums(is.na(.)) <= threshold)
  #apply function over each row
  return(exp_data_clean)
}

exp_data_filtered <- handle_na_values(exp_data)
#exp_data_filtered[, 1: 23] are control group [, 24:94] are sepsis group]

replace_na_with_median <- function (row){
  control_median <- median(row[1:23], na.rm = TRUE)
  sepsis_median <- median(row[24:94], na.rm = TRUE)
  row[1:23][is.na(row[1:23])] <- control_median
  row[24:94][is.na(row[24:94])] <- sepsis_median
  return(row)
}

exp_data_imputed <- as.data.frame(t(apply(exp_data_filtered, 1, replace_na_with_median)))

#### calculate the logfold change ####
exp_data_imputed$log_fold_change <- apply(exp_data_imputed, 1, function (row){
  control_mean <- mean(row[1:23])
  sepsis_mean <- mean(row[24:94])
  log2(sepsis_mean / control_mean)
})

diff_exp_genes <- exp_data_imputed%>%
  filter(!is.na(log_fold_change) & abs(log_fold_change) >= 2.00) %>%
  arrange(desc(log_fold_change))


write.csv(diff_exp_genes, "./Data/Diff_expr_genes", row.names = FALSE)
#sepsis vs control 
sepsis_upreg <- diff_exp_genes%>%
  dplyr::select(log_fold_change)%>%
  filter(log_fold_change > 0)

sepsis_downReg <- diff_exp_genes %>%
  dplyr::select(log_fold_change)%>%
  filter(log_fold_change < 0)

#function for GO enrichment Analysis for upregulated genes 

get_go_enrichment <- function(expr_data, GO_term){
  geneList <- expr_data$log_fold_change
  names(geneList) <- rownames(expr_data)
  #map the gene names to EntrezID
  gene_symbols <- names(geneList)
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols,
                       column = "ENTREZID", keytype = "SYMBOL",
                       multiVals = "first")
  #remove NAs and update geneList 
  valid_entrez_id <- entrez_ids[!is.na(entrez_ids)]
  geneList <- geneList[names(geneList) %in% names(valid_entrez_id)]
  names(geneList) <- valid_entrez_id
  go_result <- get_go_result(geneList, GO_term)
  return(go_result)
}

#function to get go_enrichment 
get_go_result <- function(geneList, term){
  result <- enrichGO(gene = names(geneList),
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = term,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  return(result)
}
barplot(get_go_enrichment(sepsis_upreg, "CC"), title = "GO enrichment for CC in Sepsis upregulated genes")
barplot(get_go_enrichment(sepsis_downReg, "BP"), title = "GO enrichment for BP in Sepsis Down-regulated genes")
barplot(get_go_enrichment(sepsis_downReg, "CC"), title = "GO enrichment for CC in Sepsis Down-regulated genes")


#There is no significant enrichment of Molecular Functions and Biological Processes
#for the upregulated genes in sepsis 
# and there are no significant enrichment for Molecular Functions and Cellular Components for 
#down-regulated sepsis genes 

#now let's do same thing for differential gene expression in control 
exp_data_control <- exp_data_imputed%>%
  dplyr::select(-log_fold_change)
exp_data_control$log_fold_change <- apply(exp_data_control, 1, function(row){
  control_mean <- mean(row[1:23])
  sepsis_mean <- mean(row[24:94])
  log2(control_mean/sepsis_mean)
})
#now filtered upregulated or downregulated genes 
diff_exp_control_genes <- exp_data_control%>%
  filter(!is.na(log_fold_change) & abs(log_fold_change) >= 2.00) %>%
  arrange(desc(log_fold_change))

#head(diff_exp_control_genes)
control_upreg <- diff_exp_control_genes%>%
  dplyr::select(log_fold_change)%>%
  filter(log_fold_change > 0)
control_downReg <- diff_exp_control_genes%>%
  dplyr::select(log_fold_change)%>%
  filter(log_fold_change < 0)

#Now barplots for GO Enrichment in control genes 
barplot(get_go_enrichment(control_upreg, "BP"), title = "GO enrichment for BP in Control upregulated genes")
barplot(get_go_enrichment(control_upreg, "CC"), title = "GO enrichment for CC in Control upregulated genes")

barplot(get_go_enrichment(control_downReg, "MF"), title = "GO enrichment for MF in Control down_regulated genes")
barplot(get_go_enrichment(control_downReg, "CC"), title = "GO enrichment for CC in Control down_regulated genes")

