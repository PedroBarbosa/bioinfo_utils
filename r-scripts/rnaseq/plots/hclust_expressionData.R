#####################################################
################Feature counts read counts ##########
#####################################################
args = commandArgs(trailingOnly=TRUE)
#features_counts_data=read.table(args[1], header=TRUE,sep="\t")
features_counts_data=read.table("~/Downloads/christian/genes_feature_counts.txt", header=TRUE,sep="\t")

#Remove gene_name col
df_numbers = features_counts_data[ , -which(names(features_counts_data) %in% c("gene_name"))]

#Set gene id as rowname
rownames(df_numbers) = df_numbers[,'Geneid']
df_numbers[,'Geneid'] <- NULL
colnames(df_numbers) <- sub('Aligned.sortedByCoord.out.bam', '', colnames(df_numbers))
#df_numberic = mutate_all(df_numbers, function(x) as.numeric(as.character(x)))

#Log 2 transformation
log_matrix<-log(df_numbers, base = 2)
log_matrix[mapply(is.infinite, log_matrix)] <- 0.00000000000001

pdf("clustering_rawReads.pdf")
clusters <- hclust(dist(t(df_numbers)))
#clusters <- hclust(dist(t(log_matrix)))
plot(clusters)
dev.off()

#K-means trial - dirty
clusters <- kmeans(df_numbers, 5)
plot(t(df_numbers), col = clusters$cluster)
points(clusters$centers,col = 1:3, pch = 8, cex=2)

#####################################################
##################STRINGTIE GENE COUNTS #############
#####################################################
list_files <- list.files(pattern="*_1.txt")

colClasses_FPKM <- c(NA, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", NA, "NULL")
colClasses_TPM <- c(NA, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", NA)
df_tpm <- lapply(list_files, read.table, sep = "\t", colClasses = colClasses_TPM) 
expr_measure <- "TPM"

#Initialize final dataframe from the first sample in the list
final_df = df_tpm[[1]]
#Rename columns
names(final_df) <- as.character(unlist(final_df[1,]))
final_df <- final_df[-1,]
names(final_df)[names(final_df) == expr_measure] <- gsub('.txt', '', list_files[[1]])

#Remove gene that in some samples appears twice
#n_occur <- data.frame(table(df$`Gene ID`))
#print(n_occur[n_occur$Freq > 1,])
final_df <- final_df[- grep("ENSG00000275395.5", final_df$`Gene ID`),]

#Iterate over list of dfs to merge TPM from all the samples
i = 1
for (df in df_tpm) {
  if (i != 1){
    print(paste("Processing ", list_files[[i]]))
    #Rename cols
    names(df) <- as.character(unlist(df[1,]))
    df <- df[-1,]
    names(df)[names(df) == expr_measure] <- gsub('.txt', '', list_files[[i]])
  
    #Remove gene that in some samples appears twice
    df <- df[- grep("ENSG00000275395.5", df$`Gene ID`),]
    
    #Merge dfs
    final_df=merge(final_df, df, by="Gene ID")
   
  }
  i = i + 1
}
#Set gene ID as rownames
rownames(final_df) <- df[,1]
final_df[,1] <- NULL

pdf("clustering_TPM.pdf")
clusters <- hclust(dist(t(final_df)))
plot(clusters)
dev.off()

#################################################
############PAIRWISE CORRELATIONS################
#################################################
library(corrplot)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(gridExtra)
colnames(final_df) <- lapply(colnames(final_df), function(x) substr(x,1,nchar(x)-2))
final_df <- mutate_all(final_df, function(x) as.numeric(as.character(x)))
corrplot(cor(final_df))
write.table(format(cor(final_df,method = "pearson"), digits = 4), file = "all_vs_all.csv", dec = ".", sep = "\t")
pairwise_comparisons <- vector(mode="list", length=9)
pairwise_comparisons[[1]] = c("GD1", "GD2")
pairwise_comparisons[[2]] = c("GD1", "GD3")
pairwise_comparisons[[3]] = c("GD2", "GD3")
pairwise_comparisons[[4]] = c("S_3y_CBE", "S_10y_CBE")
pairwise_comparisons[[5]] = c("S_3y_CBE", "S_11y_CBE")
pairwise_comparisons[[6]] = c("S_10y_CBE", "S_11y_CBE")
pairwise_comparisons[[7]] = c("S_3y_CTRL", "S_10y_CTRL")
pairwise_comparisons[[8]] = c("S_3y_CTRL", "S_11y_CTRL")
pairwise_comparisons[[9]] = c("S_10y_CTRL", "S_11y_CTRL")

##Multiple plot
plotlist=list()
i=1
for (pair in pairwise_comparisons) {
    p<-ggplot(final_df, aes(x=final_df[,pair[1]], y=final_df[, pair[2]])) +
      geom_point(color='brown') +
      labs(x=pair[1], y=pair[2]) +
      stat_cor(method="pearson")
    
    plotlist[[i]] = ggplotGrob(p)
    i= i + 1 
}

grid.arrange(grobs=plotlist,nrow=6,ncol=3,bottom="Pairwise comparisons of gene expression values (TPMs)", 
             gp=gpar(fontface="bold", col="YELLOW", fontsize=15))
