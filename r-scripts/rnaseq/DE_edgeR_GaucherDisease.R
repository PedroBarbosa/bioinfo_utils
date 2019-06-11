library(edgeR)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(RColorBrewer)

#Input from nextflow
args = commandArgs(trailingOnly=TRUE)
#features_counts_data=read.table(args[1], header=TRUE,sep="\t")
features_counts_data=read.table("~/Downloads/christian/featureCounts/genes_feature_counts.txt", header=TRUE,sep="\t")

#n_occur <- data.frame(table(features_counts_data$`gene_name`))
#print(n_occur[n_occur$Freq > 1,])

#Remove gene_name col
df_numbers = features_counts_data[ , -which(names(features_counts_data) %in% c("gene_name"))]

#Set gene id as rowname and reorder columns
rownames(df_numbers) = df_numbers[,'Geneid']
df_numbers[,'Geneid'] <- NULL
df_numbers <- df_numbers[c(6,7,1,2,8,3,4,5,9)]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

###############################################################
############################edgeR func #######################
edgeR_func <- function(data, group, group_order) {
  d0 <- DGEList(counts = data, group=group)
  d0 <- calcNormFactors(d0)
  d0 <- estimateCommonDisp(d0,verbose=TRUE)
  d0 <- estimateTagwiseDisp(d0)
  de.tgw <- exactTest(d0, pair=group_order)
  summary(decideTestsDGE(de.tgw, p.value=0.01))
  data$twd <- d0$tagwise.dispersion
  hist(data$twd, breaks=20, xlim=c(0,3))
  data <- cbind(data, de.tgw$table)
  data$FDR <- p.adjust(method="fdr",p=data$PValue)
  final <- data[data$FDR<0.01,]
  final <- final[order(final$logFC, decreasing = TRUE),, drop=FALSE]
  final <- rownames_to_column(final, "Geneid")
  data <- rownames_to_column(data,'Geneid')
  return(list(final, data))
}

###############################################################
##############PATIENTS (GD) versus HEALTHY CONTROLS############
###############################################################
#Remove CBE samples and very low expressed genes
df <- df_numbers[, -grep("CBE", colnames(df_numbers))]
dim(df)

#More than half of the listed genes are not expressed
df_subset <- df[rowSums(df)>10,]
dim(df_subset)

#Set group
group <- c(rep("Patients",3),rep("Control",3))

#Call the edgeR function
list_out <- edgeR_func(df_subset, group, c("Control","Patients"))
final=list_out[[1]]
final$gene_name <- features_counts_data$gene_name[match(final$Geneid, features_counts_data$Geneid)]
genedesc <- getBM(attributes=c('external_gene_name','description'), 
                  filters = 'external_gene_name', values = final$gene_name, mart =ensembl)
final$description <- genedesc$description[match(final$gene_name,genedesc$external_gene_name)]
write.table(final[,c("Geneid","gene_name","description","logFC","PValue","FDR")], quote= FALSE, sep="\t",file="GD_vs_Controls.csv",row.names=FALSE)
DE_genes_GD_CTRL <- final$gene_name
DE_genes_GD_CTRL_UP <- final[final$logFC > 0,]$gene_name
DE_genes_GD_CTRL_DOWN <- final[final$logFC < 0,]$gene_name

all_data=list_out[[2]]
all_data$gene_name <- features_counts_data$gene_name[match(all_data$Geneid, features_counts_data$Geneid)]
zeb2_GD_vs_controls =dplyr::filter(all_data, grepl("GBA",gene_name))

########################################################################
#######################CBE vs Healthy controls##########################
########################################################################

####################edgeR################################3
df <- df_numbers[c(1,2,3,7,8,9,4,5,6)]
df <- df[, -grep("GD", colnames(df))]
dim(df)

#More than half of the listed genes are not expressed
df_subset <- df[rowSums(df)>10,]
dim(df_subset)

#Set group
group <- c(rep("CBE_treatment",3),rep("Control",3))

#Call edgeR
list_out <- edgeR_func(df_subset, group, c("Control", "CBE_treatment"))
final_CBE_control_edgeR=list_out[[1]]
final_CBE_control_edgeR$gene_name <- features_counts_data$gene_name[match(final_CBE_control_edgeR$Geneid, features_counts_data$Geneid)]
genedesc <- getBM(attributes=c('external_gene_name','description'), 
                  filters = 'external_gene_name', values = final_CBE_control_edgeR$gene_name, mart =ensembl)
final_CBE_control_edgeR$description <- genedesc$description[match(final_CBE_control_edgeR$gene_name,genedesc$external_gene_name)]
write.table(final_CBE_control_edgeR[,c("Geneid","gene_name","description","logFC","PValue","FDR")], quote= FALSE, sep="\t",file="CBE_vs_Controls.csv",row.names=FALSE)
DE_genes_CBE_CTRL <- final_CBE_control_edgeR$gene_name
DE_genes_CBE_CTRL_UP <-  final_CBE_control_edgeR[final_CBE_control_edgeR$logFC > 0,]$gene_name
DE_genes_CBE_CTRL_DOWN <-  final_CBE_control_edgeR[final_CBE_control_edgeR$logFC < 0,]$gene_name
DE_genes_CBE_CTRL_DF <- final_CBE_control_edgeR[,c("Geneid", "gene_name")]

all_data=list_out[[2]]
all_data$gene_name <- features_counts_data$gene_name[match(all_data$Geneid, features_counts_data$Geneid)]
zeb2_CBE_vs_controls =dplyr::filter(all_data, grepl("GBA",gene_name))


#############DESEQ2###############
setwd("~/Downloads/christian/featureCounts/")
source("~/git_repos/bioinfo_utils/r-scripts/rnaseq/exploratory_rna_seq.R")
features_counts_data=read.table("genes_feature_counts.txt", header=TRUE,sep="\t")
df_numbers = features_counts_data[ , -which(names(features_counts_data) %in% c("gene_name"))]

#Set gene id as rowname and reorder columns
rownames(df_numbers) = df_numbers[,'Geneid']
df_numbers[,'Geneid'] <- NULL

df_numbers <- df_numbers[, c("GD1", "GD2", "GD3", "S_3y_CTRL", "S_10y_CTRL", "S_11y_CTRL", "S_3y_CBE", "S_10y_CBE", "S_11y_CBE")]
groups <- as.factor(rep(c("GD","CTRL", "CBE"), each=3))
coldata <- data.frame(groups, row.names=colnames(df_numbers))
dds <- DESeqDataSetFromMatrix(countData = df_numbers,
                              colData = coldata,
                              design = ~ groups)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$groups<-relevel(dds$groups,ref="CTRL")
log2cutoff <- 1
padjcutoff <- 0.01

group_combination <- c("groups","CBE", "CTRL")
list_de <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38")
all_annot_cbe <- annotate_results(list_de[[1]], "hg38") 
write.table(all_annot_cbe, quote= FALSE, row.names = TRUE,  sep="\t",file="CBE_vs_Controls_all_annotated.csv")

##FGSEA###
source("~/git_repos/bioinfo_utils/r-scripts/rnaseq/exploratory_rna_seq.R")
fgseaResTidy <- run_fgsea_analysis(all_annot_cbe, "~/Downloads/h.all.v6.2.symbols.gmt", "Hallmark_pathways")
fgseaResTidy <- run_fgsea_analysis(all_annot_cbe, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")


##GO over representation analysis
goseq.results <- run_goseq_analysis(all_annot_cbe, "hg19", rownames(list_de[[2]]))
head(goseq.results[goseq.results$term == "negative chemotaxis",])

upset(fromList(list(GD_vs_CTR_edgeR=final_CBE_control_edgeR$gene_name, GD_vs_CTR_deseq2=final_CBE_control_deseq2$gene_name)),
      empty.intersections = "on", text.scale = 2)#, order.by = "freq")

plotMA(resLFC, ylim=c(-2,2), xlim=c(0,10))
write.table(final[,c("gene_name","description","log2FoldChange","pvalue","padj")], quote= FALSE, sep="\t",file="CBE_vs_Controls_Deseq2.csv")

######################################################################
###########################GD vs CBE treatment########################
######################################################################
df <- df_numbers
df <- df[, -grep("CTRL", colnames(df))]
dim(df)

#More than half of the listed genes are not expressed
df_subset <- df[rowSums(df)>10,]
dim(df_subset)

#Set group
group <- c(rep("Patients",3),rep("CBE_treatment",3))

#Call edgeR
list_out <- edgeR_func(df_subset, group, c("CBE_treatment", "Patients"))
final=list_out[[1]]
final$gene_name <- features_counts_data$gene_name[match(final$Geneid, features_counts_data$Geneid)]
genedesc <- getBM(attributes=c('external_gene_name','description'), 
                  filters = 'external_gene_name', values = final$gene_name, mart =ensembl)
final$description <- genedesc$description[match(final$gene_name,genedesc$external_gene_name)]
write.table(final[,c("Geneid","gene_name","description","logFC","PValue","FDR")], quote= FALSE, sep="\t",file="GD_vs_CBE.csv",row.names=FALSE)
DE_genes_GD_CBE <- final$gene_name
DE_genes_GD_CBE_UP <- final[final$logFC >0,]$gene_name
DE_genes_GD_CBE_DOWN <- final[final$logFC < 0,]$gene_name

all_data=list_out[[2]]

all_data$gene_name <- features_counts_data$gene_name[match(all_data$Geneid, features_counts_data$Geneid)]



################INTERSECTION#################
library(UpSetR)

length(DE_genes_GD_CTRL_UP)
length(DE_genes_GD_CTRL_DOWN)
length(DE_genes_CBE_CTRL_UP)
length(DE_genes_CBE_CTRL_DOWN)
length(DE_genes_GD_CBE_UP)
length(DE_genes_GD_CBE_DOWN)

upset(fromList(list(GD_vs_CTRL=DE_genes_GD_CTRL, CBE_vs_CTRL=DE_genes_CBE_CTRL,
                    GD_vs_CBE=DE_genes_GD_CBE)), empty.intersections = "on", text.scale = 2)#, order.by = "freq")

upset(fromList(list(GD_vs_CTRL_UP = DE_genes_GD_CTRL_UP, GD_vs_CTRL_DOWN = DE_genes_GD_CTRL_DOWN,
                    CBE_vs_CTRL_UP = DE_genes_CBE_CTRL_UP, CBE_vs_CTRL_DOWN = DE_genes_CBE_CTRL_DOWN,
                    GD_vs_CBE_UP = DE_genes_GD_CBE_UP, GD_vs_CBE_DOWN = DE_genes_GD_CBE_DOWN)), text.scale = 2)#, order.by = "freq")

GD_CTRL_vs_GD_CBE <- intersect(DE_genes_GD_CTRL, DE_genes_GD_CBE)
GD_CTRL_vs_CBE_CTRL <- intersect(DE_genes_GD_CTRL, DE_genes_CBE_CTRL)
GD_CBE_vs_CBE_CTRL <- intersect(DE_genes_GD_CBE, DE_genes_CBE_CTRL)

intersect(DE_genes_GD_CTRL_DOWN, DE_genes_CBE_CTRL_UP)
intersect(DE_genes_GD_CTRL_UP, DE_genes_CBE_CTRL_DOWN)
genes_GD_CTRL_intercept_GD_CBE <- getBM(attributes=c('external_gene_name','description'), 
                  filters = 'external_gene_name', values = GD_CTRL_vs_GD_CBE, mart =ensembl)

genes_GD_CTRL_intercept_CBE_CTRL <- getBM(attributes=c('external_gene_name','description'), 
                                        filters = 'external_gene_name', values = GD_CTRL_vs_CBE_CTRL, mart =ensembl)

genes_GD_CBE_intercept_CBE_CTRL <- getBM(attributes=c('external_gene_name','description'), 
                                        filters = 'external_gene_name', values = GD_CBE_vs_CBE_CTRL, mart =ensembl)
setwd("~/Downloads/christian/")
write.table(genes_GD_CTRL_intercept_GD_CBE, quote= FALSE, sep="\t",file="GD_CTRL_intercept_GD_CBE",row.names=FALSE)
write.table(genes_GD_CTRL_intercept_CBE_CTRL, quote= FALSE, sep="\t",file="GD_CTRL_intercept_CBE_CTRL",row.names=FALSE)
write.table(genes_GD_CBE_intercept_CBE_CTRL, quote= FALSE, sep="\t",file="GD_CBE_intercept_CBE_CTRL",row.names=FALSE)


###############################################################################
#################### TPMs INPUT TO PRODUCE HEATMAP ############################
###############################################################################
setwd("/home/pbarbosa/Downloads/christian/stringtie_TPM")
#list gene expression counts from stringtie files
list_files <- list.files(pattern="*gene_abund.txt")

#cols to import: since it is RF rna-seq, will extract the reverse read counts per file
cols_TPM <- c(NA,NA, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", NA)
df_counts <- lapply(list_files, read.table, sep = "\t", colClasses = cols_TPM) 

#Initialize final dataframe from the first sample in the list
final_df = df_counts[[1]]
final_df = final_df[-1,]
#Rename columns
names(final_df) <- c("geneID", "geneName",gsub('_1Aligned.sortedByCoord.out.gene_abund.txt','', list_files[[1]]))

#Remove gene that in some samples appears twice
n_occur <- data.frame(table(final_df$geneID))
print(n_occur[n_occur$Freq > 1,])
final_df <- final_df[!final_df$geneID == 'ENSG00000275395.5',]

#Iterate over list of dfs to merge counts from all the samples
i = 1
for (df in df_counts) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Remove gene name col and rename the remainings
    df[2] <- NULL
    names(df) <- c("geneID",gsub('_1Aligned.sortedByCoord.out.gene_abund.txt','', list_files[[i]]))
    
    #Remove dups
    df <- df[!df$geneID == 'ENSG00000275395.5',]
    #Merge dfs
    final_df=merge(final_df, df, by="geneID")
  }
  i = i + 1
}

####### Total distinct genes with DE in all comparisons ######
length(as.vector.factor(DE_genes_GD_CTRL))
gd_ctrl <- as.vector.factor(DE_genes_GD_CTRL)

length(as.vector.factor(DE_genes_CBE_CTRL))
cbe_ctrl <- as.vector.factor(DE_genes_CBE_CTRL)

length(as.vector.factor(DE_genes_GD_CBE))
gd_cbe <- as.vector.factor(DE_genes_GD_CBE)

ALL_UNIQUE_DE_GENES <- c(gd_ctrl, cbe_ctrl, gd_cbe)
ALL_UNIQUE_DE_GENES <- unique(ALL_UNIQUE_DE_GENES)
length(ALL_UNIQUE_DE_GENES)
ALL_UNIQUE_DE_GENES <- cbe_ctrl

###### Filter TPM matrix by those genes to produce heatmap #####
final_df <- final_df %>% remove_rownames %>% column_to_rownames(var="geneID")
filtered_df <- final_df %>%
  rowwise() %>% 
  filter(any(ALL_UNIQUE_DE_GENES %in% geneName))

filtered_df <- filtered_df %>% remove_rownames %>% column_to_rownames(var="geneName")
genenames <- rownames(filtered_df)
filtered_df <- mutate_all(filtered_df, function(x) as.numeric(as.character(x)))

filtered_df <- log(filtered_df,2)
filtered_df[mapply(is.infinite, filtered_df)] <- 0.00000000000001
filtered_df <- as.data.frame(filtered_df)
rownames(filtered_df) <- genenames


##### Plot heatmap #######
library(heatmap3)
colors=brewer.pal(n = 8, name = "Blues")
plotheatmap <- function(x) {
  heatmap3(x,
            #main = "Genes correlation", # heat map title
            #density.info="none",  # turns off density plot inside color legend
            #trace="none",         # turns off trace lines inside the heat map
            #margins =c(3,5),     # widens margins around plot
            col=colors,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            #dendrogram="col",     # only draw a row dendrogram
            Colv=FALSE, 
            Rowv = FALSE, # turn off column clustering
            showRowDendro = FALSE,
            scale="row",
            cexRow=0.1,
            cexCol=0.8,
            breaks = c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3)
  )
}
pdf("TPM_heatmap_DE_genes.pdf")
plotheatmap(filtered_df)
write.table(filtered_df, quote= FALSE, sep="\t", file="DE_genes_log2_tpm.txt",row.names=TRUE, col.names = TRUE)
dev.off()



###############################################################
#######################GO and GSEA analysis ###################
###############################################################
setwd("~/Downloads/christian/GO_analysis/")

##########CBE vs CONTROLS###########
referenceGeneList <- rownames(df_numbers)
referenceGeneList <- gsub("\\..*", "", referenceGeneList)

#####from edgeR#####
toGO_analysis_all_edgeR <- final_CBE_control_edgeR$Geneid
toGO_analysis_all_edgeR <- gsub("\\..*", "", toGO_analysis_all_edgeR)
write.table(toGO_analysis_all_edgeR, quote= FALSE, file="CBE_CTRL_genes_all_edgeR.txt",row.names=FALSE, col.names = FALSE)

toGO_analysis_up_edgeR <- final_CBE_control_edgeR[final_CBE_control_edgeR$logFC > 0,]$Geneid
toGO_analysis_up_edgeR <- gsub("\\..*", "", toGO_analysis_up_edgeR)
write.table(toGO_analysis_up_edgeR, quote= FALSE, file="CBE_CTRL_genes_up_edgeR.txt",row.names=FALSE, col.names = FALSE)

toGO_analysis_down_edgeR <- final_CBE_control_edgeR[final_CBE_control_edgeR$logFC < 0,]$Geneid
toGO_analysis_down_edgeR <- gsub("\\..*", "", toGO_analysis_down_edgeR)
write.table(toGO_analysis_down_edgeR, quote= FALSE, file="CBE_CTRL_genes_down_edgeR.txt",row.names=FALSE, col.names = FALSE)

write.table(referenceGeneList, quote= FALSE, file="referenceGeneList_ensembl.txt", col.names = FALSE,row.names=FALSE)


####from DESeq2####
toGO_analysis_all_deseq2 <- rownames(final_CBE_control_deseq2)
toGO_analysis_all_deseq2 <- gsub("\\..*", "", toGO_analysis_all_deseq2)
write.table(toGO_analysis_all_deseq2, quote= FALSE, file="CBE_CTRL_genes_all_deseq2.txt",row.names=FALSE, col.names = FALSE)

toGO_analysis_up_deseq2 <- rownames(subset(final_CBE_control_deseq2, log2FoldChange > 0))
toGO_analysis_up_deseq2 <- gsub("\\..*", "", toGO_analysis_up_deseq2)
write.table(toGO_analysis_up_deseq2, quote= FALSE, file="CBE_CTRL_genes_up_deseq2.txt",row.names=FALSE, col.names = FALSE)

toGO_analysis_down_deseq2 <- rownames(subset(final_CBE_control_deseq2, log2FoldChange < 0))
toGO_analysis_down_deseq2 <- gsub("\\..*", "", toGO_analysis_down_deseq2)
write.table(toGO_analysis_down_deseq2, quote= FALSE, file="CBE_CTRL_genes_down_deseq2.txt",row.names=FALSE, col.names = FALSE)

write.table(referenceGeneList, quote= FALSE, file="referenceGeneList_ensembl.txt", col.names = FALSE,row.names=FALSE)


