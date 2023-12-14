#### load the packages ####
library(tidyverse)
library(ggplot2)
library(GENIE3)
library(readr)
library(dplyr)
library(tibble)
library(doParallel)
library(doRNG)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#### load data ####
data <- read_csv('./Data/GDS2808.soft.csv', show_col_types = FALSE)
human_tfs <- read_csv('./Data/human_tfs.csv', show_col_types = FALSE)

exp_data <- data[!duplicated(data$IDENTIFIER),]
exp_data <- exp_data%>%
  dplyr::select(-'ID_REF')%>%
  column_to_rownames(var = 'IDENTIFIER')

#exp_data <- na.omit(data) #--> 148 rows 
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
# expression values are the fluorescent signal converted to digital signal
# that provides the measure of gene-expression level.

replace_na_with_median <- function (row){
  control_median <- median(row[1:23], na.rm = TRUE)
  sepsis_median <- median(row[24:94], na.rm = TRUE)
  row[1:23][is.na(row[1:23])] <- control_median
  row[24:94][is.na(row[24:94])] <- sepsis_median
  return(row)
}

exp_data_imputed <- as.data.frame(t(apply(exp_data_filtered, 1, replace_na_with_median)))

exp_data_imputed$log_fold_change <- apply(exp_data_imputed, 1, function (row){
  control_mean <- mean(row[1:23])
  sepsis_mean <- mean(row[24:94])
  log2(sepsis_mean / control_mean)
})

diff_exp_genes <- exp_data_imputed%>%
  filter(!is.na(log_fold_change) & abs(log_fold_change) >= 2.00) %>%
  arrange(desc(log_fold_change))

diff_exp_genes <- diff_exp_genes%>%
  dplyr::select(-log_fold_change)
#separate the control and sepsis samples 
control_samples <- as.matrix(diff_exp_genes[, 1:23])
sepsis_samples <- as.matrix(diff_exp_genes[, 24:94])
label_data <- read_csv('./Data/labels.csv', col_names = TRUE, show_col_types = FALSE)


#function to get regulators for GENIE3
get_regulators <- function(data, tfs){
  target_genes <- row.names(data)
  is_regulator <- target_genes %in% tfs$Symbol
  regulatoryGenes <- target_genes[is_regulator]
}

control_regulatoryGens <- get_regulators(control_samples, human_tfs)
sepsis_regulatoryGenes <- get_regulators(sepsis_samples, human_tfs)

#apply genie3 for gene regulatory analysis 
#function to run genie3

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

control_linklist <- get_linkList(control_samples, control_regulatoryGens)
control_toplink <- head(control_linklist, 50)
network_graph_control <- graph_from_data_frame(control_toplink, directed = TRUE)
plot(network_graph_control, vertex.size = 5, 
     vertex.label.cex = 0.8,, edge.arrow.size = 0.5,  layout = layout_with_fr)

sepsis_linklist <- get_linkList(sepsis_samples, sepsis_regulatoryGenes)
sepsis_toplink <- head(sepsis_linklist, 50)
network_graph_sepsis <- graph_from_data_frame(sepsis_toplink, directed = TRUE)
plot(network_graph_sepsis, vertex.size = 5, 
     vertex.label.cex = 0.8,, edge.arrow.size = 0.5,  layout = layout_with_fr)

