library(ggplot2)
library("stringr")
library(ggbeeswarm)
library(beeswarm)
library(shinycssloaders)
libs <- c("tidyverse", "ggVennDiagram", "BiocManager",
          "DESeq2", "edgeR", "limma","gplots")

for (package in libs) {
  suppressPackageStartupMessages(require(package, 
                                         quietly = T, 
                                         character.only = T))
  require(package, character.only = T)
}

deseq_normalize <- function(count_data, meta_data, design_formula) {
  deseq <- DESeqDataSetFromMatrix(
    countData=select(count_data, -c(gene_id)),
    colData=meta_data,
    design=~1
  )
  deseq <- estimateSizeFactors(deseq)
  result <- counts(deseq, normalized=TRUE)
  result <- as_tibble(result) %>%
    mutate(gene_id=count_data$gene_id) %>%
    relocate(gene_id)
  return(result)
}
meta_info_from_labels <- function(sample_names) {
  timepoint <- sapply(sample_names, timepoint_from_sample)
  treatment <- sapply(sample_names, sample_treatment)
  result <- tibble(
    sample=sample_names,
    timepoint=timepoint,
    treatment=treatment
  )
  return(result)
}

#getting the counts data
data <- readr::read_csv("GSE150450_gene_count_matrix.csv.gz")

metadata1 <- meta_info_from_labels(colnames(data)[colnames(data) != "gene_id"])

deseq_data <- deseq_normalize(data, metadata1, ~1)
#write.csv(deseq_data,"deseq_data.csv")