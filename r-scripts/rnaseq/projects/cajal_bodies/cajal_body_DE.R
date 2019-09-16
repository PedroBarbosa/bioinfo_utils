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
setwd("~/Downloads/cajal_body_project/")
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","gene_name","Chr","Start","End","Strand","Length")]
colnames(read_counts) <- sub('X.home.pedro.barbosa.cb_data.star_alignments.trim_galored.','', colnames(read_counts))
colnames(read_counts) <- sub('_Aligned.sortedByCoord.out.bam','', colnames(read_counts))
colnames(read_counts) <- sub('myC_time','t', colnames(read_counts))

#Filter genes where at least five reads are present in at least 2 samples
filter <- apply(read_counts, 1, function(x) length(x[x>5])>=2)
filtered <- read_counts[filter,]



###################################################
########  RUVseq spike ins control   ##############
###################################################
timepoints <- as.factor(rep(c("0", "12", "24", "48"), each=3))
coldata <- data.frame(timepoints, row.names=colnames(filtered))
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = coldata,
                              design = ~ timepoints)

genes <- rownames(dds)[grep("^ENS", rownames(dds))]
spikes <- rownames(dds)[grep("^ERCC", rownames(dds))]
set <- newSeqExpressionSet(counts(dds))

#Data exploration
colors <- brewer.pal(8, "Set2")
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(set, col=colors, cex=1.2)

#Remove variation using spike-in control genes
#Estimate factors of unwanted variation
normalized <- RUVg(set, spikes, k=1)
plot_umwanted_variation_factors(pData(normalized), dds$timepoints, 1)

#Plot normalized counts
plotRLE(normCounts(normalized), outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(normCounts(normalized), col=colors, cex=1.2)

#rerun DEseq with the new design to reestimante the parametes and results
dds$W1 <- normalized$W_1
#dds$W2 <- normalized$W_2
#design(dds) <- ~ W2 + W1 + timepoints
design(dds) <- ~ W1 + timepoints

log2cutoff <- 1
padjcutoff <- 0.05

#########################################
############# 12 vs 0 ###################
#########################################
group_combination <- c("timepoints","12", "0")
dds$timepoints <-relevel(dds$timepoints, ref="0")
list_de_12_0 <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38", explore_data = FALSE)
rownames(list_de_12_0[[1]]) <- gsub("\\..*", "", rownames(list_de_12_0[[1]]))
all_annot_12_0 <- merge(as(list_de_12_0[[1]], "data.frame") , ensembl_genes, all.x = T, by=0)
#all_annot_12_0 <- annotate_results(list_de_8_0[[1]], "hg38") 

write.table(all_annot_12_0, quote= FALSE, row.names = TRUE, sep="\t",file="T12_vs_T0_all_annotated.csv")

fgseaResTidy <- run_fgsea_analysis(all_annot_12_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_12_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_12_0, "hg19", list_de_12_0[[2]]$Row.names)

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

fgseaResTidy <- run_fgsea_analysis(all_annot_48_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_48_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_48_0, "hg19", list_de_48_0[[2]]$Row.names)



########################################
#############CAJAL BODY GENES ##########
########################################
cb_genes <- read.table("~/Downloads/cajal_body_project/CB_genes.txt")
histone_genes <- read.table("~/Downloads/cajal_body_project/Histone_genes.txt")

de_genes_12_0 <- list_de_12_0[[2]]$symbol
cb_genes <- unlist(cb_genes, use.names=FALSE)
histone_genes <- unlist(histone_genes, use.names=FALSE)
paste("Number of de genes 12 vs 0:", length(de_genes_12_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_12_0)))
paste("Number of cajal body genes in the whole gene list:", length(which(cb_genes %in% all_annot_12_0$symbol)))
#cb_genes not in tables: AK6 FAM71B FBL FRG1 FRG1P LSM2 NELFE NHP2L1 PPP1R10 SMN TERT
cb_genes_present_12_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_12_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_12_0)))
paste("Number of histone genes in the whole gene list:", length(which(histone_genes %in% all_annot_12_0$symbol)))
histone_genes_present_12_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_12_0)]))
#histone genes not in tables: All there

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

df_cb_genes <- as.data.frame(t(bind_rows(lapply(list(cb_genes_present_12_0, 
                                                     cb_genes_present_24_0, cb_genes_present_48_0), as.data.frame.list))))
colnames(df_cb_genes) <- c("12_vs_0", "24_vs_0", "48_vs_0")

df_histone_genes <-  as.data.frame(t(bind_rows(lapply(list(histone_genes_present_12_0, 
                                                           histone_genes_present_24_0, histone_genes_present_48_0), as.data.frame.list))))
colnames(df_histone_genes) <- c("12_vs_0", "24_vs_0", "48_vs_0")
write.table(df_cb_genes, quote= FALSE, row.names = FALSE, sep="\t",file="CB_genes_DE.csv")
write.table(df_histone_genes, quote= FALSE, row.names = FALSE, sep="\t",file="histone_genes_DE.csv")
