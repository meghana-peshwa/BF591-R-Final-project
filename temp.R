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
#getting the counts data
data <- readr::read_csv("GSE150450_gene_count_matrix.csv.gz")
first_column <- as.list(data[, 1]$gene_id)
first_column
genes <- paste(shQuote(first_column), collapse=", ")
#getting the metadata
metadata <- read.table("SraRunTable.txt",sep=",")
names(metadata) <- metadata[1,]
metadata <- metadata[-1,]

#metadata_updated = metadata[,-1]
metadata_updated = metadata[,-5]
metadata_updated = metadata_updated[,-11]
metadata_updated = metadata_updated[,-11]
metadata_updated = metadata_updated[,-18]
metadata_updated = metadata_updated[,-18]
column_names = names(metadata_updated)
column_type = vector()
unique_values = vector()
metadata_updated <- transform(metadata_updated, AvgSpotLen = as.numeric(AvgSpotLen), 
               Bases = as.numeric(Bases),
               Bytes = as.numeric(Bytes))
for(i in names(metadata_updated)){
  column_type = append(column_type,class(metadata_updated[[i]]))
  if(class(metadata_updated[[i]])=="character") {
    unique_values <- append(unique_values,paste(unique(metadata_updated[[i]]), collapse=";"))
  }
  else{
    unique_values <- append(unique_values,mean(metadata_updated[[i]]))
  }
}
sample_info <- tibble(
  column_names=column_names,
  column_type=column_type,
  unique_values=unique_values
)

metadata_updated$lifestage <- factor(metadata_updated$lifestage, labels=c("adult", "larvae"))
metadata_updated$timepoint <- factor(metadata_updated$timepoint, labels=c("high2","high3","low2","low3","high1","low1"))
metadata_updated$treatment <- factor(metadata_updated$treatment, labels=c("control","fluctuating"))
metadata_updated$sex <- factor(metadata_updated$sex, labels=c("male","female",""))

barplot(metadata_updated$lifestage,
        main="Multiple Bar Plots",
        ylab="Count"
)

timepoint_from_sample <- function(str) {
  print(substr(str,3,4))
  return(substr(str,3,4))
}

sample_treatment <- function(str) {
  return(substr(str,1,1))
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


metadata1 <- meta_info_from_labels(colnames(data)[colnames(data) != "gene_id"])

#deseq_data <- deseq_normalize(data, metadata1, ~1)
#write.csv(deseq_data,"deseq_data.csv")
deseq_data <- read.csv("../deseq_data.csv")[,-c(1)]
counts_data_filtered <- deseq_data[rowSums(deseq_data[-1]!=0)>=80,]
data <- counts_data_filtered[,-c(1)]
variance <- apply(data, 1, var)
counts_data_filtered_final <- counts_data_filtered[variance>=quantile(variance,95/100,na.rm=TRUE),]

filter_counts <- function(counts_data, zero_filter, variance_filter) {
  if(is.null(counts_data))
    return(NULL)
  counts_data_filtered <- counts_data[rowSums(counts_data[-1]!=0)>=zero_filter,]
  data <- counts_data_filtered[,-c(1)]
  variance <- apply(data, 1, var)
  counts_data_filtered_final <- counts_data_filtered[variance>=quantile(variance,variance_filter/100,na.rm=TRUE),]
  summary <- tibble(
    no_of_samples=length(counts_data_filtered)-1,
    total_number_of_genes=nrow(counts_data),
    percent_genes_passing_filter=nrow(counts_data_filtered_final)/nrow(counts_data)*100,
    percent_genes_not_passing_filter=(nrow(counts_data)-nrow(counts_data_filtered_final))/nrow(counts_data)*100,
  )
  return(list(counts_data_filtered_final,summary))
}

counts_data_filtered_final <- filter_counts(deseq_data,80,95)

plot_variance_vs_median <- function(data, scale_y_axis=FALSE, title="") {
  median <- apply(data[, -c(ncol(data))], 1, median)
  variance <- apply(data[, -c(ncol(data))], 1, var)
  plot_data <- tibble(median=median, variance=variance, filtered = data$filtered)
  plot_data$rank <- rank(plot_data$median)
  plot <- ggplot(plot_data, aes(x=rank, y=variance)) +
    geom_point(mapping=aes(x=rank, y=variance, color=filtered)) +
    xlab("Rank(median)") +
    ylab("Variance") +
    ggtitle(title)
  if (scale_y_axis) {
    plot <- plot + scale_y_log10()
  }
  return(plot)
}

df1 <- deseq_data %>%
  left_join(counts_data_filtered_final %>% transmute(gene_id, filtered = 'no')) %>%
  replace_na(list(filtered = 'yes'))

df2 <- df1[,-c(1)]
plot_variance_vs_median(df2,TRUE,"")

plot_zeros_vs_median <- function(data, title="") {
  zeros <- rowSums(data[, -c(ncol(data))] == 0)
  median <- apply(data[, -c(ncol(data))], 1, median)
  plot_data <- tibble(zeros=zeros, median=median, filtered = data$filtered)
  plot_data$rank <- rank(plot_data$median)
  plot <- ggplot(plot_data, aes(x=rank, y=zeros)) +
    geom_point(mapping=aes(x=rank, y=zeros, color=filtered)) +
    xlab("Rank(median)") +
    ylab("Number of zeros") + 
    ggtitle(title)
  return(plot)
}

plot_zeros_vs_median(df2,"")
zeros <- rowSums(df2[, -c(ncol(data))] == 0)
median <- apply(df2[, -c(ncol(df2))], 1, median)

matrix_heatmap <- as.matrix(counts_data_filtered_final[[1]][,-c(1)])
heatmap.2(log10(matrix_heatmap), trace="none", key=TRUE, scale="row")


plot_pca <- function(data, meta, pc1, pc2, title="") {
  pca <- prcomp(t(data))
  meta$PC1 <- pca$x[ , pc1]
  meta$PC2 <- pca$x[ , pc2]
  variance <- pca$sdev^2 / sum( pca$sdev^2 )
  pca_plot <- ggplot(meta, aes(x=PC1, y=PC2, col=lifestage)) +
    geom_point() +
    xlab(paste0("PC",toString(pc1),": ",round(variance[1] * 100),"% variance")) +
    ylab(paste0("PC",toString(pc2),": ",round(variance[2] * 100),"% variance")) +
    ggtitle(title)
  return(pca_plot)
}

plot_pca(deseq_data[c(-1)],result,2,4,"")

run_deseq <- function(count_dataframe, coldata, count_filter) {
  filter_counts<-as.matrix(count_dataframe[rowSums(count_dataframe)>=count_filter,])
  dds <- DESeqDataSetFromMatrix(countData = round(filter_counts),
                                colData = coldata,
                                design = formula(~ treatment+lifestage))
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}
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

draw_table <- function(dataf, slider) {
  if(is.null(dataf)) {
    return(NULL)
  }
  print(dataf)
  print(slider)
  dataf <- dataf[dataf$padj<(1*10^slider),]
  print('hi')
  #dataf$pvalue <- formatC(dataf$pvalue,digits = 10)
  #dataf$padj <- formatC(dataf$padj,digits = 10)
  print(dataf)
  return(dataf)
}
table <- draw_table(deseq2_results,-150)
dataf<-na.omit(deseq2_results)
dataf <- dataf[dataf$padj<(10^(-15)),]
dataf<-na.omit(dataf)
if(nrow(dataf)==0) {
  print('heyyy')
}

deseq_data <- read.csv("deseq_data.csv",header=TRUE)
deseq_data <- deseq_data[,-c(1)]
gene_data <- t(deseq_data[deseq_data$gene_id == 'FBgn0036214',][,-c(1)])

timepoint_from_sample <- function(str) {
  return(substr(str,3,4))
}

sample_treatment <- function(str) {
  return(substr(str,1,1))
}

sex_from_sample <- function(str) {
  return(substr(str,length(str),length(str)))
}

lifestage_from_sample <- function(str) {
  return(substr(str,2,2))
}

cols <- colnames(deseq_data)[colnames(deseq_data) != "gene_id"]
timepoint <- sapply(cols, timepoint_from_sample)
treatment <- sapply(cols, sample_treatment)
sex <- sapply(cols, sex_from_sample)
lifestage <- sapply(cols, lifestage_from_sample)
result <- as.data.frame(tibble(
  sample=cols,
  timepoint=c(timepoint),
  treatment=c(treatment),
  sex=c(sex),
  lifestage=c(lifestage)
))
x = "timepoint"
barplot(table(result[,x]),
        main=paste0("Bar chart for ",x),
        ylab="Count")
beeswarm(counts ~ timepoint,data=result,method="swarm",pch = 16, col = rainbow(25))
ggplot(result, aes(x=timepoint,y=counts,color=timepoint)) +
  geom_beeswarm()
v = "treatment"
plot <- ggplot(result, aes(x=result[,v],y=counts,fill=result[,v]))+xlab(v) +
  ylab("Counts")
plot+geom_boxplot()

library(ggridges)

ggplot(metadata_updated) +
  geom_violin(aes(x=treatment,y=Bases,fill=treatment))

