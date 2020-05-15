library(DESeq2)
library(tidyverse)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

###############################
#######  FROM SALMON  #########
###############################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/christian/human_fibroblasts/npa_npc_april_2020/salmon/")
files <- list.files(path=".", pattern = "*sf")

groups <- as.factor(unlist(lapply(strsplit(files, split = "_"), `[`, 2)))
samples <-  as.factor(unlist(lapply(strsplit(files, split = "_"), `[`, 1)))
names(files) <- names(files) <- samples
samples_df <- data.frame(groups, samples)

tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38.txt")
tx_ids <- readr::read_tsv(file = "13y_CTRL_GM02037_quant.sf") %>% .$Name 
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID version`, `Gene stable ID`) , by = c("value" = "Transcript stable ID version")) 

txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene)


############################
###### DESeq2 dataset ######
############################
coldata = data.frame(timepoints, row.names =  colnames(txi_salmon$counts))
dds_salmon <- DESeqDataSetFromTximport(txi_salmon,
                                       colData = coldata,
                                       design = ~ timepoints)

keep <- rowSums(counts(dds_salmon)) >= 10
dds_salmon <- dds_salmon[keep,]