library(dplyr)
library(tidyr)
library(reshape2)
library(corrplot)
library(tidyverse)
library(tibble)
library(ggplot2)
library(ggpubr)
library(gridExtra)
setwd("~/Downloads/cajal_body_project/")

##################################################################
######################SPIKE IN CHECK##############################
##################################################################
####MIX######
spike_in_table=read.table("spike_in_concentrations.txt", sep="\t",header = TRUE)
mix1=subset(spike_in_table,select=c("ERCC.ID","concentration.in.Mix.1..attomoles.ul."))
mix1[3] = log(mix1[2],2)
colnames(mix1)=c("ERCC","Mix1_concentration","Log2_concentration")
mix1[, 2:3] <- sapply(mix1[, 2:3], as.numeric)
#mix_final <- mix1 %>% remove_rownames %>% column_to_rownames(var="ERCC")
#mix <- mutate_all(mix1[], function(x) as.numeric(as.character(x)))


#################################################################
################FROM STAR RAW COUNTS#############################
#################################################################
ERCC_star_counts=read.table("reverse_ERCC_star_counts.txt")
ERCC_star_counts_split <-ERCC_star_counts %>% separate(V1, c("sample", "ERCC"), ":")
colnames(ERCC_star_counts_split)[3] <- "#Unstranded"
colnames(ERCC_star_counts_split)[4] <- "#Forward"
colnames(ERCC_star_counts_split)[5] <- "#Reverse"

reshaped=data.frame(dcast(ERCC_star_counts_split, ERCC ~ sample , value.var = '#Reverse'))
df_final <- mutate_all(reshaped[,-1], function(x) as.numeric(as.character(x)))
df_final <- log(df_final,2)
df_final[mapply(is.infinite, df_final)] <- 0.00000001
df_final['ERCC'] = reshaped['ERCC']
MERGED=merge(x = mix1, y = df_final, by = "ERCC", all = TRUE)
MERGED <-MERGED %>% remove_rownames %>% column_to_rownames(var="ERCC")
corrplot(cor(MERGED[,-c(1,2)], method = c("pearson")))


#################################################################
###################FROM FEATURE COUNTS###########################
#################################################################
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","gene_name","Chr","Start","End","Strand","Length")]
colnames(read_counts) <- sub('X.home.pedro.barbosa.cb_data.star_alignments.trim_galored.','', colnames(read_counts))
colnames(read_counts) <- sub('_Aligned.sortedByCoord.out.bam','', colnames(read_counts))
colnames(read_counts) <- sub('myC_time','t', colnames(read_counts))
log_read_counts <- log(read_counts, 2)
log_read_counts[mapply(is.infinite, log_read_counts)] <- 0.00000001
spikes <- log_read_counts[grep("^ERCC", rownames(log_read_counts)),]
mix1 <-mix1 %>% remove_rownames %>% column_to_rownames(var="ERCC")
MERGED=merge(x = mix1, y = spikes, by = 0, all = TRUE)
MERGED <-MERGED %>% remove_rownames %>% column_to_rownames(var="Row.names")

##Multiple plot
plotlist=list()
i=1
for (col in colnames(MERGED)) {

  if (!col == "Log2_concentration" && !col == "Mix1_concentration"){
    print(col)
    p<-ggplot(MERGED, aes(x=MERGED$Log2_concentration, y=MERGED[,col])) +
            geom_point(color='brown') +
            labs(y=col) +
            theme(axis.title.x=element_blank()) +
            stat_cor(method="spearman")

    plotlist[[i]] = ggplotGrob(p)
    i= i + 1 
  }

}
grid.arrange(grobs=plotlist,nrow=4,ncol=3,bottom="Log2(ERCC concentration)", left="Log2(read counts)",gp=gpar(fontface="bold", col="YELLOW", fontsize=15))

#Per sample spike-in correlation
write.table(format(cor(MERGED[,-c(1,2)],method = "pearson"), digits = 4), file = "all_vs_all.csv", dec = ".", sep = "\t")



######################################################################
####################Diff Expr Analysis Workflow#######################
######################################################################
library(RUVSeq)
library(DESeq2)
library("apeglm")
library("pheatmap")
library("RColorBrewer")
library(biomaRt)
################################################################
##############FROM STAR COUNTS WITH 1 FILE PER SAMPLE###########
################################################################
setwd("~/Downloads/cajal_body_project/star_counts/")

#list gene counts from STAR file
list_files <- list.files(pattern="*ReadsPerGene.out.tab")
list_files
#cols to import: since it is RF rna-seq, will extract the reverse read counts per file
cols <- c(NA,"NULL", "NULL", NA)
df_counts <- lapply(list_files, read.table, sep = "\t", colClasses = cols) 

#Initialize final dataframe from the first sample in the list
final_df = df_counts[[1]]
final_df = final_df[-c(1:4),]
#Rename columns
names(final_df) <- c("geneID",gsub('_ReadsPerGene.out.tab','', list_files[[1]]))

#Remove gene that in some samples appears twice
#n_occur <- data.frame(table(final_df$geneID))
#print(n_occur[n_occur$Freq > 1,])

#Iterate over list of dfs to merge counts from all the samples
i = 1
for (df in df_counts) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Rename cols
    colnames(df) <- c("geneID",gsub('_ReadsPerGene.out.tab','', list_files[[i]]))

    #Merge dfs
    final_df=merge(final_df, df, by="geneID")
  }
  i = i + 1
}

#Fix rownames and data types
rownames(final_df) = final_df[,'geneID']
final_df[,'geneID'] <- NULL

########################################################
##########FROM FEATURE COUNTS ##########################
########################################################
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","gene_name","Chr","Start","End","Strand","Length")]
colnames(read_counts) <- sub('X.home.pedro.barbosa.cb_data.star_alignments.trim_galored.','', colnames(read_counts))
colnames(read_counts) <- sub('_Aligned.sortedByCoord.out.bam','', colnames(read_counts))
colnames(read_counts) <- sub('myC_time','t', colnames(read_counts))

#Filter genes where at least five reads are present in at least 2 samples
filter <- apply(read_counts, 1, function(x) length(x[x>5])>=2)
filtered <- read_counts[filter,]


#####################################################################
################RNA seq data exploration using DESeq2 ###############
#####################################################################
#The normalization steps don't take gene size into account, since it doesn't matter. 
#You do not have a bias toward longer genes, rather you have increased power to find changes in them given a constant expression level.
#This is a good thing, you do not want to try to get rid of it.

#Nice RNA seq tutorial: browseVignettes("rnaseqGene")
###############################################################
#############variance transformation func #####################
###############################################################

explore_data_based_on_transformed_variance <- function(dds, cond){
  #Deseq2 offers 2 transformation methods to stabilize variance across the mean
  #Which method to choose, good discussion in the vignette
  #This is good explore data by  clustering samples by their distance of making a PCA analysis without the bias in the variance for high counts genes
  #v_transformed <- vst(dds, blind = FALSE)
  v_transformed <- rlog(dds, blind = FALSE) 
  #v_transformed <- log2(counts(dds) + 1)

  #Sample distances
  sampleDists <- dist(t(assay(v_transformed)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(colData(v_transformed)[,cond])
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  
  #PCA
  DESeq2::plotPCA(v_transformed, intgroup = cond)
}


visualize_differences_in_variance_transformation <- function(log_t, vsd, rld) {

  #Show the effect of transformation plotting two samples against each other
  df <- bind_rows(
    as_data_frame(log_t[, 1:2]+1) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
  ggplot(df, aes(x = df$myC_time0_rep1, y = df$myC_time12_rep1)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 

}


##################################
######Differential testing########
##################################
#Recommended the DESeq function applied to raw counts,  which also takes into account the dependence of the variance
#of counts on the mean value during the dispersion estimation step.
run_DEseq_tests <- function(dds, test_contrast, coefficient, log2cutoff, padjcutoff) {
  dds_2 <- DESeq(dds)
  res <- results(dds_2, contrast=test_contrast, alpha=0.05, independentFiltering = TRUE, pAdjustMethod="BH")
  #res_24_0 <- results(dds, contrast=c("timepoints","24","0"), alpha=0.05, pAdjustMethod="BH")
  show(resultNames(res))
  if (length(resultsNames(res)) > 1) {
    res <- lfcShrink(dds_2, type="apeglm", res= res)
        
  } else {
    res <- lfcShrink(dds_2, coef = coefficient, type="apeglm", res= res)
  }
  
  res_sign <- subset(res, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
  return(list(res,res_sign))
 
}

sort_byFC <- function(res){
  sort_down <- res[ order(res$log2FoldChange), ]
  sort_up <- res[ order(res$log2FoldChange, decreasing = TRUE), ]
  return(list(sort_up,sort_down))
}


plot_single_gene_variation <- function(dds, group, gene){
  #Normalized data
  plotCounts(dds, gene = gene, intgroup=group)
  
  #From counts
  library("ggbeeswarm")
  geneCounts <- plotCounts(dds, gene = gene, intgroup = group,
                           returnData = TRUE)
  ggplot(geneCounts, aes(x = timepoints, y = count, color = group)) +
    scale_y_log10() +  geom_beeswarm(cex = 3)
}


plot_ma <- function(res_shrinked, padjcutoff){
  plotMA(res_shrinked, alpha=padjcutoff, xlab="mean of normalized counts", ylim = c(-5, 5))
  
  #with individual labels
  #topGene <- rownames(ma_12_0)[which.min(ma_12_0$padj)]
  #with(ma_12_0[topGene, ], {
  #  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  #  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
  #})
}

cluster_sign_genes <- function(dds, sigGenes, group){
  print(paste("Number of DE genes across all combinations: ", length(sigGenes)))
  print("Clustering them based on the variance transformed raw counts.")
  v_transformed <- vst(dds, blind = FALSE)
 # print(typeof)
  #counts_justSigGenes <- counts(dds)[sigGenes,]
  counts_justSigGenes <- v_transformed[sigGenes,]
  DESeq2::plotPCA(counts_justSigGenes, intgroup = group)
}

annotate_results <- function(res){
  library(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genedesc <-  getBM(attributes=c('ensembl_gene_id_version','external_gene_name','description'), 
                     filters = 'ensembl_gene_id_version', values = rownames(res), mart =ensembl)
  
  res$symbol <- genedesc$external_gene_name[match(rownames(res),genedesc$ensembl_gene_id_version)]
  res$description <- genedesc$description[match(rownames(res),genedesc$ensembl_gene_id_version)]
  return(res)
}

#RuvSEQ specific
plot_umwanted_variation_factors <-function(pdata_normalized, groups, ksize) {
  par(mfrow = c(ksize, 1), mar = c(3,5,3,1))
  for (i in 1:ksize) {
    stripchart(pdata_normalized[, i] ~ groups, vertical = TRUE, main = paste0("W", i))
    abline(h = 0)
  }
  par(mfrow = c(1, 1))
}



####################################
######Analysis ####################
####################################

#Create DDS and experimental design
timepoints <- as.factor(rep(c("0", "12", "24", "48"), each=3))
coldata <- data.frame(timepoints, row.names=colnames(filtered))
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = coldata,
                              design = ~ timepoints)

#Data exploration
explore_data_based_on_transformed_variance(dds, c("timepoints"))


#Diff expression
log2cutoff <- 1
padjcutoff <- 0.05

dds$timepoints<-relevel(dds$timepoints,ref="0")
out_12_0 <-run_DEseq_tests(dds, c("timepoints", "12", "0"), "timepoints_12_vs_0", log2cutoff, padjcutoff)
out_12_0[[2]]$symbol <- mapIds(org.Hs.eg.db,
                     keys=rownames(out_12_0[[2]]),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
rownames(out_12_0[[2]]) <- gsub("\\..*", "", rownames(out_12_0[[2]]))
out_12_0_sign_annot <- annotate_results(out_12_0[[2]])
print(paste0("Number of genes in the 12_vs_0 comparison: ", nrow(out_12_0_sign_annot)))
sorted_12_0 <- sort_byFC(res = out_12_0_sign_annot)
plot_ma(out_12_0[[1]] , padjcutoff)
#Plot one specific gene
topGene <- rownames(out_12_0[[1]])[which.min(out_12_0[[1]]$padj)]
plot_single_gene_variation(dds, group = "timepoints", topGene )

out_24_0 <-run_DEseq_tests(dds, c("timepoints", "24", "0"), "timepoints_24_vs_0", log2cutoff, padjcutoff)
out_24_0_sign_annot <- annotate_results(out_24_0[[2]])
print(paste0("Number of genes in the 24_vs_0 comparison: ", nrow(out_24_0_sign_annot)))
sorted_24_0 <- sort_byFC(res = out_24_0_sign_annot)
plot_ma(out_24_0[[1]],  padjcutoff)

out_48_0 <-run_DEseq_tests(dds, c("timepoints", "48", "0"), "timepoints_48_vs_0", log2cutoff, padjcutoff)
out_48_0_sign_annot <- annotate_results(out_48_0[[2]])
head(out_48_0_sign_annot,1)
print(paste0("Number of genes in the 48_vs_0 comparison: ", nrow(out_48_0_sign_annot)))
sorted_48_0 <- sort_byFC(res = out_48_0_sign_annot)
plot_ma(out_48_0[[1]], padjcutoff)


dds$timepoints<-relevel(dds$timepoints,ref="24")
out_48_24 <-run_DEseq_tests(dds, c("timepoints", "48", "24"), "timepoints_48_vs_24", log2cutoff, padjcutoff)
out_48_24_sign_annot <- annotate_results(out_48_24[[2]])
head(out_48_24_sign_annot,1)
sorted_48_24 <- sort_byFC(res = out_48_24_sign_annot)
print(paste0("Number of genes in the 48_vs_24 comparison: ", nrow(out_48_24_sign_annot)))
plot_ma(out_48_24[[1]], padjcutoff)

sigGenes <- unique(c(rownames(out_12_0_sign_annot), rownames(out_24_0_sign_annot), rownames(out_48_0_sign_annot)))
cluster_sign_genes(dds, sigGenes, "timepoints")

library(UpSetR)
upset(fromList(list(t12_vs_0=out_12_0_sign_annot$symbol, t24_vs_0=out_24_0_sign_annot$symbol, t48_vs_0=out_48_0_sign_annot$symbol)),
      empty.intersections = "on", text.scale = 2)#, order.by = "freq")


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

#From origianl
#set <- newSeqExpressionSet(as.matrix(filtered),
#                           phenoData = data.frame(timepoints, row.names=colnames(filtered)))
set <- newSeqExpressionSet(counts(dds))

#Data exploration
colors <- brewer.pal(8, "Set2")
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(set, col=colors, cex=1.2)

#Remove variation using spike-in control genes
#Estimate factors of unwanted variation
normalized <- RUVg(set, spikes, k=2)
plot_umwanted_variation_factors(pData(normalized), dds$timepoints, 2)

#Plot normalized counts
plotRLE(normCounts(normalized), outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(normCounts(normalized), col=colors, cex=1.2)

#rerun DEseq with the new design to reestimante the parametes and results
dds$W1 <- normalized$W_1
dds$W2 <- normalized$W_2
design(dds) <- ~ W2 + W1 + timepoints

#Data exploration
explore_data_based_on_transformed_variance(dds, c("timepoints"))

log2cutoff <- 2
padjcutoff <- 0.05
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","Chr","Start","End","Strand","Length")]

dds$timepoints<-relevel(dds$timepoints,ref="0")
out_12_0 <-run_DEseq_tests(dds, c("timepoints", "12", "0"), "timepoints_12_vs_0", log2cutoff, padjcutoff)
out_12_0_sign_annot <- out_12_0[[2]]
out_12_0_sign_annot$gene_name <- read_counts$gene_name[match(rownames(out_12_0_sign_annot), rownames(read_counts))]
#out_12_0_sign_annot <- annotate_results(out_12_0[[2]])
print(paste0("Number of genes in the 12_vs_0 comparison: ", nrow(out_12_0_sign_annot)))
sorted_12_0 <- sort_byFC(res = out_12_0_sign_annot)
plot_ma(out_12_0[[1]] , padjcutoff)
#Plot one specific gene
topGene <- rownames(out_12_0[[1]])[which.min(out_12_0[[1]]$padj)]
plot_single_gene_variation(dds, group = "timepoints", topGene )
write.table(out_12_0_sign_annot, quote= FALSE, sep="\t", file="DE_t12_vs_t0_deseq2.txt",row.names=TRUE, col.names = TRUE)

out_24_0 <-run_DEseq_tests(dds, c("timepoints", "24", "0"), "timepoints_24_vs_0", log2cutoff, padjcutoff)
out_24_0_sign_annot <- out_24_0[[2]]
out_24_0_sign_annot$gene_name <- read_counts$gene_name[match(rownames(out_24_0_sign_annot), rownames(read_counts))]
#out_24_0_sign_annot <- annotate_results(out_24_0[[2]])
print(paste0("Number of genes in the 24_vs_0 comparison: ", nrow(out_24_0_sign_annot)))
sorted_24_0 <- sort_byFC(res = out_24_0_sign_annot)
plot_ma(out_24_0[[1]],  padjcutoff)
write.table(out_24_0_sign_annot, quote= FALSE, sep="\t", file="DE_t24_vs_t0_deseq2.txt",row.names=TRUE, col.names = TRUE)

out_48_0 <-run_DEseq_tests(dds, c("timepoints", "48", "0"), "timepoints_48_vs_0", log2cutoff, padjcutoff)
out_48_0_sign_annot <- out_48_0[[2]]
out_48_0_sign_annot$gene_name <- read_counts$gene_name[match(rownames(out_48_0_sign_annot), rownames(read_counts))]
#out_48_0_sign_annot <- annotate_results(out_48_0[[2]])
print(paste0("Number of genes in the 48_vs_0 comparison: ", nrow(out_48_0_sign_annot)))
sorted_48_0 <- sort_byFC(res = out_48_0_sign_annot)
plot_ma(out_48_0[[1]], padjcutoff)
write.table(out_48_0_sign_annot, quote= FALSE, sep="\t", file="DE_t48_vs_t0_deseq2.txt",row.names=TRUE, col.names = TRUE)

dds$timepoints<-relevel(dds$timepoints,ref="24")
out_48_24 <-run_DEseq_tests(dds, c("timepoints", "48", "24"), "timepoints_48_vs_24", log2cutoff, padjcutoff)
out_48_24_sign_annot <- annotate_results(out_48_24[[2]])
head(out_48_24_sign_annot,1)
sorted_48_24 <- sort_byFC(res = out_48_24_sign_annot)
print(paste0("Number of genes in the 48_vs_24 comparison: ", nrow(out_48_24_sign_annot)))
plot_ma(out_48_24[[1]], padjcutoff)


sigGenes <- unique(c(rownames(out_12_0_sign_annot), rownames(out_24_0_sign_annot), rownames(out_48_0_sign_annot)))
cluster_sign_genes(dds, sigGenes, "timepoints")

library(UpSetR)
upset(fromList(list(t12_vs_0=rownames(out_12_0_sign_annot), t24_vs_0=rownames(out_24_0_sign_annot), t48_vs_0=rownames(out_48_0_sign_annot))),
      empty.intersections = "on", text.scale = 2)#, order.by = "freq")
upset(fromList(list(t12_vs_0=out_12_0_sign_annot$symbol, t24_vs_0=out_24_0_sign_annot$symbol, t48_vs_0=out_48_0_sign_annot$symbol)),
      empty.intersections = "on", text.scale = 2)#, order.by = "freq")


###############################################
###############Time course#####################
##############################################
dds <- DESeq(dds, test="LRT")

d