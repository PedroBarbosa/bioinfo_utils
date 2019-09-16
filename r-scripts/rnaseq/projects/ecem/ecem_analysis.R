library(tidyverse)
library(biomaRt)
#####################################################
##################SALMOM TPM COUNTS #############
#####################################################
setwd("~/Downloads/ecem/salmon/")
list_files <- list.files(pattern="*.sf")
colClasses_TPM <- c(NA, "NULL", "NULL", NA, "NULL")
df_tpm <- lapply(list_files, read.table, sep = "\t", colClasses = colClasses_TPM) 

#Initialize final dataframe from the first sample in the list
final_df = df_tpm[[1]]
#Rename columns
names(final_df) <- as.character(unlist(final_df[1,]))
final_df <- final_df[-1,]
names(final_df)[2] <- gsub('_1.fq.gz_quant.sf', '', list_files[[1]])

#Remove gene that in some samples appears twice
#n_occur <- data.frame(table(final_df$Name))
#print(n_occur[n_occur$Freq > 1,])

#Iterate over list of dfs to merge TPM from all the samples
i = 1
for (df in df_tpm) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Rename cols
    names(df) <- as.character(unlist(df[1,]))
    df <- df[-1,]
    names(df)[2] <- gsub('_1.fq.gz_quant.sf', '', list_files[[i]])
    
    #Merge dfs
    final_df=merge(final_df, df, by="Name")
    
  }
  i = i + 1
}
#Set gene ID as rownames
rownames(final_df) <- df[,1]
final_df[,1] <- NULL
final_df <- mutate_all(final_df, function(x) as.numeric(as.character(x)))

#log transform
log_final_df <- log(final_df,2)
log_final_df[mapply(is.infinite, log_final_df)] <- 0.00000001

#Clustering
#pdf("clustering_TPM_salmon_spearman_log.pdf")
#hc <- hclust(as.dist(1-cor(log_final_df, method="spearman")), method="complete") # Colums by spearman
#plot(hc)
#dev.off()

pdf("clustering_log2_TPM_salmon.pdf")
clusters <- hclust(dist(t(log_final_df)))
plot(clusters)
dev.off()



###########################
####Kallisto TPM counts ###
###########################
setwd("~/Downloads/ecem/kallisto/")
list_files <- list.files(pattern="*.tsv")
colClasses_TPM <- c(NA, "NULL", "NULL", "NULL", NA)
df_tpm <- lapply(list_files, read.table, sep = "\t", colClasses = colClasses_TPM) 

#Initialize final dataframe from the first sample in the list
final_df = df_tpm[[1]]
#Rename columns
names(final_df) <- as.character(unlist(final_df[1,]))
final_df <- final_df[-1,]
names(final_df)[2] <- gsub('_1.fq.gz_abundance.tsv', '', list_files[[1]])
#Get just trancript id
final_df[,1] <-  gsub("\\|.*", "", final_df[,1])

#Remove gene that in some samples appears twice
#n_occur <- data.frame(table(final_df$target_id))
#print(n_occur[n_occur$Freq > 1,])

#Iterate over list of dfs to merge TPM from all the samples
i = 1
for (df in df_tpm) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Rename cols
    names(df) <- as.character(unlist(df[1,]))
    df <- df[-1,]
    df[,1] <-  gsub("\\|.*", "", df[,1])
    names(df)[2] <- gsub('_1.fq.gz_abundance.tsv', '', list_files[[i]])
    #Get just trancript id
    df[,1] <-  gsub("\\|.*", "", df[,1])
    #Merge dfs
    final_df=merge(final_df, df, by="target_id")
    
  }
  i = i + 1
}

#Set gene ID as rownames
rownames(final_df) <- df[,1]
final_df[,1] <- NULL
final_df <- mutate_all(final_df, function(x) as.numeric(as.character(x)))

#log transform
log_final_df <- log(final_df,2)
log_final_df[mapply(is.infinite, log_final_df)] <- 0.00000001

#Clustering
#pdf("clustering_TPM_kallisto_spearman_log.pdf")
#hc <- hclust(as.dist(1-cor(log_final_df, method="spearman")), method="complete") # Colums by spearman
#plot(hc)
#dev.off()

pdf("clustering_log2_TPM_kallisto.pdf")
clusters <- hclust(dist(t(log_final_df)))
plot(clusters)
dev.off()


####################################
####Stringtie TPM counts############
####################################
setwd("/home/pbarbosa/Downloads/ecem/stringtie/")
#list gene expression counts from stringtie files
list_files <- list.files(pattern="*tsv")

#cols to import: since it is RF rna-seq, will extract the reverse read counts per file
cols_TPM <- c(NA,"NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", NA)
df_counts <- lapply(list_files, read.table, sep = "\t", colClasses = cols_TPM) 

#Initialize final dataframe from the first sample in the list
final_df = df_counts[[1]]
names(final_df) <- as.character(unlist(final_df[1,]))
final_df = final_df[-1,]
#Rename columns
names(final_df) <- c("geneID",gsub('_quantification.tsv','', list_files[[1]]))

#Remove gene that in some samples appears twice
n_occur <- data.frame(table(final_df$geneID))
print(n_occur[n_occur$Freq > 1,])
repeated_genes <- levels(droplevels(n_occur[n_occur$Freq > 1,]$Var1))
final_df <- final_df[!final_df$geneID %in% repeated_genes,]

#Iterate over list of dfs to merge counts from all the samples
i = 1
for (df in df_counts) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Remove gene name col and rename the remainings
    names(df) <- c("geneID",gsub('_quantification.tsv','', list_files[[i]]))
    
    #Remove dups
    n_occur <- data.frame(table(df$geneID))
    repeated_genes <- levels(droplevels(n_occur[n_occur$Freq > 1,]$Var1))
    df <- df[!df$geneID %in% repeated_genes,]
    
    #Merge dfs
    final_df=merge(final_df, df, by="geneID")
  }
  i = i + 1
}

#Set gene ID as rownames
rownames(final_df)
rownames(final_df) <- final_df[,1]
final_df[,1] <- NULL

#log transform
final_df <- mutate_all(final_df, function(x) as.numeric(as.character(x)))
log_final_df <- log(final_df,2)
log_final_df[mapply(is.infinite, log_final_df)] <- 0.00000001

#Clustering
#pdf("clustering_TPM_stringtie_spearman_log.pdf")
#hc <- hclust(as.dist(1-cor(log_final_df, method="spearman")), method="complete") # Colums by spearman
#plot(hc)
#dev.off()

pdf("clustering_log2_TPM_stringtie.pdf")
clusters <- hclust(dist(t(log_final_df), method = "euclidean"))
plot(clusters)
dev.off()



################################################
#############Feature counts for DE##############
################################################
setwd("~/Downloads/ecem/feature_counts")
features_counts_data=read.table("merged_gene_counts.txt", header=TRUE,sep="\t")

#n_occur <- data.frame(table(features_counts_data$`Geneid`))
#print(n_occur[n_occur$Freq > 1,])

#Remove gene_name col
df_numbers = features_counts_data[ , -which(names(features_counts_data) %in% c("gene_name"))]
names(df_numbers) <- c(gsub('_1Aligned.sortedByCoord.out.bam','', names(df_numbers)))
#Set gene id as rowname and reorder columns
rownames(df_numbers) = df_numbers[,'Geneid']
df_numbers[,'Geneid'] <- NULL

########################
#####Comparisons########
########################
source("~/git_repos/bioinfo_utils/r-scripts/rnaseq/exploratory_rna_seq.R")


#Filter genes where at least five reads are present in at least 2 samples
filter <- apply(df_numbers, 1, function(x) length(x[x>5])>=2)
filtered <- df_numbers[filter,]
filtered <- filtered[,c("WT_1D1", "WT_3D1", "WT_UT1", "NPC_UT1", "NPC_UT2", "NPC_1D1", "NPC_1D2", "NPC_3D1", "NPC_3D2")]

groups <- as.factor(rep(c("WT","NPC_UT", "NPC1", "NPC3"), times=c(3,2,2,2)))
coldata <- data.frame(groups, row.names=colnames(filtered))
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = coldata,
                              design = ~ groups)

log2cutoff <- 1.5
padjcutoff <- 0.01


######################################################
#####WT UT + WT 1D1 + WT3D1 vs NPC UT1 + NPC UT2######
######################################################
setwd("~/Downloads/ecem/feature_counts")
group_combination <- c("groups","NPC_UT", "WT")
DE_NPC_UT_vs_WT <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "mm10")
all_annot <- annotate_results(DE_NPC_UT_vs_WT[[1]], "mm10") 
write.table(all_annot, quote= FALSE, row.names = TRUE,  sep="\t",file="NPC_UT_vs_WT_all_annotated.csv")
all_annot <- read.table("NPC_UT_vs_WT_all_annotated.csv", row.names = 1, sep="\t", header=TRUE)
source("~/git_repos/bioinfo_utils/r-scripts/rnaseq/exploratory_rna_seq.R")
fgseaResTidy <- run_fgsea_analysis(all_annot, "~/Downloads/MousePath_All_gmt-Format.gmt", "MousePath_Hallmarks")
fgseaResTidy
fgseaResTidy <- run_fgsea_analysis(all_annot, "~/Downloads/MousePath_GO_gmt.gmt", "MousePath_GO")
##GO over representation analysis
goseq.results <- run_goseq_analysis(all_annot, "mm9", rownames(DE_NPC_UT_vs_WT[[2]]))
head(goseq.results)

####################################################
######NPC UT1 + NPC UT2 vs NPC 1D1 + NPC 1D2########
####################################################
dds$groups<-relevel(dds$groups,ref="NPC_UT")
group_combination <- c("groups","NPC1", "NPC_UT")
DE_NPC1_vs_NPC_UT <- run_analysis(dds, group_combination, 1, padjcutoff)


#########################################
#NPC UT1 + NPC UT2 vs NPC 3D1 + NPC 3D2
#########################################
group_combination <- c("groups","NPC3", "NPC_UT")
DE_NPC3_vs_NPC_UT <- run_analysis(dds, group_combination, 1, padjcutoff)


#########################################
#NPC 1D1+NPC 1D2 vs NPC 3D1+NPC 3D2######
#########################################
dds$groups<-relevel(dds$groups,ref="NPC1")
group_combination <- c("groups","NPC3", "NPC1")
DE_NPC3_vs_NPC1 <- run_analysis(dds, group_combination, 1, padjcutoff)


########################################
#WT 1D1 + WT3D1 vs NPC 1D1+NPC 1D2######
########################################
dds2 <- dds[,-3] #Remove WT_UT sample
dds2$groups<-relevel(dds2$groups,ref="WT")
group_combination <- c("groups","NPC1", "WT")
DE_NPC1_vs_WT_1_3 <- run_analysis(dds2, group_combination, log2cutoff , padjcutoff)


########################################
#WT 1D1 + WT3D1 vs NPC 3D1+NPC 3D2######
########################################
group_combination <- c("groups","NPC3", "WT")
DE_NPC3_vs_WT_1_3 <- run_analysis(dds2, group_combination, log2cutoff, padjcutoff)



l <- list(NPC1_vs_NPC_UT = DE_NPC1_vs_NPC_UT$symbol,
          NPC3_vs_NPC_UT = DE_NPC3_vs_NPC_UT$symbol,
          NPC3_vs_NPC1 = DE_NPC3_vs_NPC1$symbol)
        #NPC_UT_vs_WT = DE_NPC_UT_vs_WT$symbol,
        # NPC1_vs_WT_1_3 = DE_NPC1_vs_WT_1_3$symbol,
        #  NPC3_vs_WT_1_3 = DE_NPC3_vs_WT_1_3$symbol)

library(UpSetR)

upset(fromList(l), text.scale = 2, order.by = "freq")
     # empty.intersections = "off", text.scale = 2)#, order.by = "freq")
