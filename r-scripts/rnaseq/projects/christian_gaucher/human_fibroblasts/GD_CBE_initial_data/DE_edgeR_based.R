library(edgeR)
##############################
#######edgeR func ############
##############################
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

######################################
#########PATIENTS (GD) vs CONTROLS####
######################################
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