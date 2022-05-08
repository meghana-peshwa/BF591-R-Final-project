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

run_deseq <- function(count_dataframe, coldata, count_filter) {
  filter_counts<-as.matrix(count_dataframe[rowSums(count_dataframe)>=count_filter,])
  dds <- DESeqDataSetFromMatrix(countData = round(filter_counts),
                                colData = coldata,
                                design = formula(~ treatment+lifestage))
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}
deseq_data <- read.csv("./deseq_data.csv")[,-c(1)]
coldata1 <- data.frame(
  sample = colnames(deseq_data)[colnames(deseq_data) != "gene_id"],
  treatment = factor(c(sapply(colnames(deseq_data)[colnames(deseq_data) != "gene_id"], sample_treatment))),
  lifestage = factor(c(sapply(colnames(deseq_data)[colnames(deseq_data) != "gene_id"], lifestage_from_sample)))
)
rownames(coldata1) <- NULL
temp <- deseq_data[,-c(1)]
deseq2_results<-run_deseq(temp,coldata1,10)
deseq2_results <- data.frame(deseq2_results)

#write.csv(deseq2_results,"deseq2_results.csv")