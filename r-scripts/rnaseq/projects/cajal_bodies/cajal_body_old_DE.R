library(DESeq2)
library(RUVSeq)
library(dplyr)
library(tidyr)
library(reshape2)
library(corrplot)
library(tidyverse)
library(tibble)
library(ggplot2)
library(ggpubr)
library(gridExtra)
source("~/git_repos/bioinfo_utils/r-scripts/rnaseq/exploratory_rna_seq.R")

########################################################
##########FROM FEATURE COUNTS ##########################
########################################################
setwd("~/Downloads/cajal_body_project_old/")
ensembl_genes <- read.table("~/Downloads/ensembl_genes_with_description.txt", quote = "", row.names = 1,sep="\t", header=TRUE)
names(ensembl_genes)[names(ensembl_genes) == "Gene.name"] <- "symbol"
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","gene_name","Chr","Start","End","Strand","Length")]
colnames(read_counts) <- sub('X.home.pedro.barbosa.cb_old_data.star.','', colnames(read_counts))
colnames(read_counts) <- sub('Aligned.sortedByCoord.out.bam','', colnames(read_counts))
read_counts <- read_counts[, c(7, 8, 2, 10, 4, 9, 1, 5, 3, 6)]

#Filter genes where at least five reads are present in at least 2 samples
filter <- apply(read_counts, 1, function(x) length(x[x>5])>=2)
filtered <- read_counts[filter,]

###################################################
########  RUVseq spike ins control   ##############
###################################################
timepoints <- as.factor(rep(c("WT", "0", "8", "24", "48"), each=2))
timepoints <-relevel(timepoints, ref="WT")
coldata <- data.frame(timepoints, row.names=colnames(filtered))
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = coldata,
                              design = ~ timepoints)


log2cutoff <- 1
padjcutoff <- 0.05

#########################################
############# 8 vs 0 ###################
#########################################
group_combination <- c("timepoints","8", "0")
dds$timepoints <-relevel(dds$timepoints, ref="0")
list_de_8_0 <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38", explore_data = FALSE)
rownames(list_de_8_0[[1]]) <- gsub("\\..*", "", rownames(list_de_8_0[[1]]))
all_annot_8_0 <- merge(as(list_de_8_0[[1]], "data.frame") , ensembl_genes, all.x = T, by=0)
#all_annot_8_0 <- annotate_results(list_de_8_0[[1]], "hg38") 


write.table(all_annot_8_0, quote= FALSE, row.names = TRUE, sep="\t",file="T8_vs_T0_all_annotated.csv")

fgseaResTidy <- run_fgsea_analysis(all_annot_8_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_8_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_8_0, "hg19", list_de_8_0[[2]]$Row.names)

#########################################
############# 24 vs 0 ###################
#########################################
group_combination <- c("timepoints","24", "0")
list_de_24_0 <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38", explore_data = FALSE)
rownames(list_de_24_0[[1]]) <- gsub("\\..*", "", rownames(list_de_24_0[[1]]))
all_annot_24_0 <- merge(as(list_de_24_0[[1]], "data.frame") , ensembl_genes, all.x = T, by=0)
#all_annot_24_0 <- annotate_results(list_de_24_0[[1]], "hg38") 
write.table(all_annot_24_0, quote= FALSE, row.names = TRUE, sep="\t",file="T24_vs_T0_all_annotated.csv")
#all_annot_24_0 <- read.table("T24_vs_T0_all_annotated.csv", row.names = 1, sep="\t", header=TRUE)

fgseaResTidy <- run_fgsea_analysis(all_annot_24_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_24_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_24_0, "hg19", list_de_24_0[[2]]$Row.names)

#########################################
############# 48 vs 0 ##################
#########################################
group_combination <- c("timepoints","48", "0")
list_de_48_0 <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38", explore_data = FALSE)
#all_annot_48_0 <- annotate_results(list_de_48_0[[1]], "hg38") 
rownames(list_de_48_0[[1]]) <- gsub("\\..*", "", rownames(list_de_48_0[[1]]))
all_annot_48_0 <- merge(as(list_de_48_0[[1]], "data.frame") , ensembl_genes, all.x = T, by=0)
write.table(all_annot_48_0, quote= FALSE, row.names = TRUE, sep="\t",file="T48_vs_T0_all_annotated.csv")
#all_annot_48_0 <- read.table("T48_vs_T0_all_annotated.csv", row.names = 1, sep="\t", header=TRUE)

fgseaResTidy <- run_fgsea_analysis(all_annot_24_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_24_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_24_0, "hg19", list_de_24_0[[2]]$Row.names)


########################################
#############CAJAL BODY GENES ##########
########################################
cb_genes <- read.table("~/Downloads/cajal_body_project/CB_genes.txt")
histone_genes <- read.table("~/Downloads/cajal_body_project/Histone_genes.txt")

de_genes_8_0 <- list_de_8_0[[2]]$symbol
cb_genes <- unlist(cb_genes, use.names=FALSE)
histone_genes <- unlist(histone_genes, use.names=FALSE)
paste("Number of de genes 8 vs 0:", length(de_genes_8_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_8_0)))
paste("Number of cajal body genes in the whole gene list:", length(which(cb_genes %in% all_annot_8_0$symbol)))
#cb_genes not in tables: FAM71B FRG1P  NHP2L1 SMN   TERT
cb_genes_present_8_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_8_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_8_0)))
paste("Number of histone genes in the whole gene list:", length(which(histone_genes %in% all_annot_8_0$symbol)))
histone_genes_present_8_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_8_0)]))
#histone genes not in tables: CPSF73 TROVE2

de_genes_24_0 <- list_de_24_0[[2]]$symbol
paste("Number of de genes 24 vs 0:", length(de_genes_24_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_24_0)))
cb_genes_present_24_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_24_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_24_0)))
histone_genes_present_24_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_24_0)]))


de_genes_48_0 <- list_de_48_0[[2]]$symbol
paste("Number of de genes 48 vs 0:", length(de_genes_48_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_48_0)))
cb_genes_present_48_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_48_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_48_0)))
histone_genes[histone_genes %in% de_genes_48_0]
histone_genes_present_48_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_48_0)]))

df_cb_genes <- as.data.frame(t(bind_rows(lapply(list(cb_genes_present_8_0, 
                                                     cb_genes_present_24_0, cb_genes_present_48_0), as.data.frame.list))))
colnames(df_cb_genes) <- c("8_vs_0", "24_vs_0", "48_vs_0")

df_histone_genes <-  as.data.frame(t(bind_rows(lapply(list(histone_genes_present_8_0, 
                                                           histone_genes_present_24_0, histone_genes_present_48_0), as.data.frame.list))))
colnames(df_histone_genes) <- c("8_vs_0", "24_vs_0", "48_vs_0")
write.table(df_cb_genes, quote= FALSE, row.names = FALSE, sep="\t",file="old_data_CB_genes_DE.csv")
write.table(df_histone_genes, quote= FALSE, row.names = FALSE, sep="\t",file="old_data_histone_genes_DE.csv")

