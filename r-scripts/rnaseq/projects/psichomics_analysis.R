library(psichomics)
library(tidyverse)
library(tidylog)
library(writexl)
library(ggplot2)
library(ggrepel)
library(limma)
setwd("/Users/pbarbosa/analysis/yvan_rna_seq/psichomics/star_counts")

fixInNamespace("prepareGeneQuantSTAR", pos="package:psichomics")
#setnames(table, colnames(table)[[2]], paste0("col", index))

###############################
#########PREPARE DATA##########
###############################
prepareGeneQuant("Ctrl_F1505921"="F1505921_CtrlReadsPerGene.out.tab",
                 "Ctrl_F1505922"="F1505922_CtrlReadsPerGene.out.tab",
                 "Ctrl_F1505923"="F1505923_CtrlReadsPerGene.out.tab",
                 "Ctrl_F1505924"="F1505924_CtrlReadsPerGene.out.tab",
                 "Ctrl_F1505925"="F1505925_CtrlReadsPerGene.out.tab",
                 "DCM_F1505926"="F1505926_DCMReadsPerGene.out.tab",
                 "ICM_F1505927"="F1505927_ICMReadsPerGene.out.tab",
                 "DCM_F1505928"="F1505928_DCMReadsPerGene.out.tab",
                 "DCM_F1505929"="F1505929_ICMReadsPerGene.out.tab",
                 "ICM_F1505930"="F1505930_ICMReadsPerGene.out.tab",
                 "ICM_F1505931"="F1505931_ICMReadsPerGene.out.tab",
                 "DCM_F1505932"="F1505932_DCMReadsPerGene.out.tab",
                 "ICM_F1505933"="F1505933_ICMReadsPerGene.out.tab",
                 "DCM_F1505934"="F1505934_DCMReadsPerGene.out.tab",
                 "ICM_F1505935"="F1505935_ICMReadsPerGene.out.tab",
                 "DCM_F1505936"="F1505936_DCMReadsPerGene.out.tab",
                 "DCM_F1505937"="F1505937_DCMReadsPerGene.out.tab",
                 "DCM_F1505938"="F1505938_DCMReadsPerGene.out.tab",
                 "ICM_F1505939"="F1505939_ICMReadsPerGene.out.tab",
                 "ICM_F1505940"="F1505940_ICMReadsPerGene.out.tab",
                 "ICM_F150594"="F1505941_ICMReadsPerGene.out.tab",
                 "ICM_F1505942"="F1505942_ICMReadsPerGene.out.tab",
                 "DCM_F1505943"="F1505943_DCMReadsPerGene.out.tab",
                 "ICM_F1505944"="F1505944_ICMReadsPerGene.out.tab",
                 "DCM_F1505945"="F1505945_DCMReadsPerGene.out.tab",
                 "DCM_F1505946"="F1505946_DCMReadsPerGene.out.tab",
                 strandedness = "stranded (reverse)")

prepareJunctionQuant("Ctrl_F1505921"="F1505921_CtrlSJ.primary_chrom.out.tab",
                     "Ctrl_F1505922"="F1505922_CtrlSJ.primary_chrom.out.tab",
                     "Ctrl_F1505923"="F1505923_CtrlSJ.primary_chrom.out.tab",
                     "Ctrl_F1505924"="F1505924_CtrlSJ.primary_chrom.out.tab",
                     "Ctrl_F1505925"="F1505925_CtrlSJ.primary_chrom.out.tab",
                     "DCM_F1505926"="F1505926_DCMSJ.primary_chrom.out.tab",
                     "ICM_F1505927"="F1505927_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505928"="F1505928_DCMSJ.primary_chrom.out.tab",
                     "DCM_F1505929"="F1505929_ICMSJ.primary_chrom.out.tab",
                     "ICM_F1505930"="F1505930_ICMSJ.primary_chrom.out.tab",
                     "ICM_F1505931"="F1505931_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505932"="F1505932_DCMSJ.primary_chrom.out.tab",
                     "ICM_F1505933"="F1505933_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505934"="F1505934_DCMSJ.primary_chrom.out.tab",
                     "ICM_F1505935"="F1505935_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505936"="F1505936_DCMSJ.primary_chrom.out.tab",
                     "DCM_F1505937"="F1505937_DCMSJ.primary_chrom.out.tab",
                     "DCM_F1505938"="F1505938_DCMSJ.primary_chrom.out.tab",
                     "ICM_F1505939"="F1505939_ICMSJ.primary_chrom.out.tab",
                     "ICM_F1505940"="F1505940_ICMSJ.primary_chrom.out.tab",
                     "ICM_F150594"="F1505941_ICMSJ.primary_chrom.out.tab",
                     "ICM_F1505942"="F1505942_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505943"="F1505943_DCMSJ.primary_chrom.out.tab",
                     "ICM_F1505944"="F1505944_ICMSJ.primary_chrom.out.tab",
                     "DCM_F1505945"="F1505945_DCMSJ.primary_chrom.out.tab",
                     "DCM_F1505946"="F1505946_DCMSJ.primary_chrom.out.tab"
)

###################
####LOAD DATA######
###################
data <- loadLocalFiles("input_data/")
setwd("input_data/")
geneExpr <- data$Data$`Gene expression`
junctionQuant <- data$Data$`Junction quantification`
groups <- list("Ctrl"=grep("Ctrl", colnames(geneExpr), value=TRUE), 
               "DCM"=grep("DCM", colnames(geneExpr), value=TRUE), 
               "ICM"=grep("ICM", colnames(geneExpr), value=TRUE))
human <- listSplicingAnnotations()[[3]]
annotation <- loadAnnotation(human)

##################
####SPLICING #####
##################
psi <- quantifySplicing(annotation, junctionQuant, minReads=10)
#####################
####Preprocessing####
#####################
psi <- psi[rowSums(is.na(psi)) <= 10,]
events <- rownames(psi)
psi <- psi %>%
  mutate(variance=rowVars(., na.rm = T)) 
rownames(psi) <- events
psi <- psi[ psi$variance > quantile(psi$variance , 0.90, na.rm = T) , ]
psi <- psi[, !(colnames(psi) %in% c("variance"))]
psi_pca <- performPCA(t(psi))
plotPCA(psi_pca, groups = groups, loadings = FALSE)

#manual pca
percentage <- round(psi_pca$sdev / sum(psi_pca$sdev) * 100, 2)
psi_pca_manual <- as.data.frame(psi_pca$x)
psi_pca_manual$group <- sapply( strsplit(as.character(row.names(psi_pca_manual)), "_"), "[[", 1 )
percentage <- paste( colnames(psi_pca_manual), "(", paste0( as.character(percentage), "%", ")", sep="") )

p <- ggplot(psi_pca_manual,aes(x=PC1, y=PC2, color=group))
p<-p+geom_point(size=3) + xlab(percentage[1]) + ylab(percentage[2])

###############
##Exploration##
###############
table <- calculateLoadingsContribution(psi_pca)
head(table, 5)
event <- "AFE_18_+_34710445_34679592_34755976_DTNA"
info  <- queryEnsemblByEvent(event, species="human", assembly="hg38")
plotTranscripts(info, event = event)
plotDistribution(psi["AFE_18_+_34710445_34679592_34755976_DTNA",], psi = T, rug = T, groups = groups)


######################
######ALL GROUPS######
######################
diffSplicing_all <- diffAnalyses(psi, groups = groups, analyses = "kruskal")
deltaPSIthreshold_all <- diffSplicing_all$`Kruskal p-value (BH adjusted)` < 0.05
delta_psi_all <- diffSplicing_all[which(deltaPSIthreshold_all), ]
dim(delta_psi_all)
#No events. Speak with nuno

#######################
######DCM vs Ctrl #####
#######################
dcm_ctrl <- groups[c("DCM", "Ctrl")]
diffSplicing_dcm_ctrl <- diffAnalyses(psi, groups = dcm_ctrl, analyses = "wilcoxRankSum")
deltaPSI_dcm_ctrl <- subset(diffSplicing_dcm_ctrl,diffSplicing_dcm_ctrl$`Wilcoxon p-value (BH adjusted)`< 0.1
                            & abs(diffSplicing_dcm_ctrl$`∆ Median`) > 0.1)  

final_dcm_ctrl <- deltaPSI_dcm_ctrl %>% rownames_to_column("event_id") %>%
  select(.,c(1, 5, 2, 9, 10, 12, 13, 20, 21, 23)) %>% 
  arrange_at(ncol(.), desc) %>%
  write_excel_csv(., "DCM_vs_CTRL.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")
write_xlsx(final_dcm_ctrl, path  = "DCM_vs_CTRL.xlsx", col_names=T, format_headers = T)


#######################
#######ICM vs Ctrl ####
#######################
icm_ctrl <- groups[c("ICM", "Ctrl")]
diffSplicing_icm_ctrl <- diffAnalyses(psi, groups = icm_ctrl, analyses = "wilcoxRankSum")
deltaPSI_icm_ctrl <- subset(diffSplicing_icm_ctrl, diffSplicing_icm_ctrl$`Wilcoxon p-value (BH adjusted)` < 0.1 & 
                              abs(diffSplicing_icm_ctrl$`∆ Median`) > 0.1)   

final_icm_ctrl <- deltaPSI_icm_ctrl %>% rownames_to_column("event_id") %>%
  select(.,c(1, 5, 2, 9, 10, 12, 13, 20, 21, 23)) %>% 
  arrange_at(ncol(.), desc) %>%
  write_excel_csv(., "ICM_vs_CTRL.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")
write_xlsx(final_icm_ctrl, path  = "ICM_vs_CTRL.xlsx", col_names=T, format_headers = T)

#######################
##CONCAT###############
#######################
final_df <- bind_rows(final_dcm_ctrl, final_icm_ctrl) %>%
  write_excel_csv(., "all_concat_vs_Ctrl.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")
write_xlsx(final_df, path  = "all_concat_vs_Ctrl.xlsx", col_names=T, format_headers = T)

#######################
####GENE EXPRESSION####
#######################
geneExprFiltered <- geneExpr[rowSums(geneExpr)>10,]
geneExprNorm <- normaliseGeneExpression(geneExprFiltered, log2transform=TRUE)
pcaGE_all <- performPCA(t(geneExprNorm))
plotPCA(pcaGE_all, groups=groups)

#Set gene ID as rownames
rownames(final_df) <- df[,1]
final_df[,1] <- NULL

clusters <- hclust(dist(t(geneExprNorm)))
plot(clusters)

###################
####manual pca#####
###################
ge_pca <- prcomp(t(geneExprNorm), center = T,scale. = T)
percentage <- round(ge_pca$sdev / sum(ge_pca$sdev) * 100, 2)
ge_pca <- as.data.frame(ge_pca$x)
ge_pca$group <- sapply( strsplit(as.character(row.names(ge_pca)), "_"), "[[", 1 )
percentage <- paste( colnames(ge_pca), "(", paste0( as.character(percentage), "%", ")", sep="") )

p <- ggplot(ge_pca,aes(x=PC1, y=PC2, color=group))
p<-p+geom_point(size=3) + xlab(percentage[1]) + ylab(percentage[2])
p

#############
####groups###
#############
to_match <- c("DCM", "ICM")
groups_sick_vs_ctrl <- list("Ctrl"=grep("Ctrl", colnames(geneExpr), value=TRUE), 
               "Patient"= grep("Ctrl", colnames(geneExpr), invert = T, value=TRUE))

pgenes <- rownames(geneExprNorm)
plotDistribution(geneExprNorm["ENSG00000163359", ], groups, psi=FALSE)


ge_sick_ctrl <- geneExprNorm[ , unlist(groups_sick_vs_ctrl), drop=FALSE]
isFromGroup1 <- colnames(ge_sick_ctrl) %in% groups_sick_vs_ctrl[[1]]
design <- cbind(1, ifelse(isFromGroup1, 0, 1))
fit <- lmFit(ge_sick_ctrl, design = design)
ebayes_fit <- eBayes(fit, trend=T)
pvalueAdj <- "BH"
summary <- topTable(ebayes_fit, number=nrow(fit), coef=2, sort.by = "none", adjust.method = pvalueAdj, confint = T)
names(summary) <- c("log2_Fold_change", "conf_int1", "conf_int2", "moderated_t-statistics", "avg_expression","p-value", 
                    paste0("p-value(",pvalueAdj,"_adjusted)"), "B-statistics")
attr(summary, "groups") <- groups_sick_vs_ctrl
stats <- diffAnalyses(ge_sick_ctrl, groups_sick_vs_ctrl, "basicStats", geneExpr = T, pvalueAdjust = NULL)
final <- cbind(stats, summary)


#Statistical significance
logFCthreshold  <- abs(final$`log2_Fold_change`) > 1
pvalueThreshold <- final$`p-value(BH_adjusted)` < 0.01

genes_match <- read_tsv("/Users/pbarbosa/analysis/genome_utilities/hg38/mart_hg38.txt")
cardio_genes <- read_tsv("/Users/pbarbosa/analysis/yvan_rna_seq/cardiomyopathy_gene_panel_ids.txt", col_names = F) %>%
  mutate(gene_id = gsub("\\..*","",.[[1]])) %>%
  select(gene_id) %>% left_join(genes_match, by=c("gene_id" = "Gene stable ID")) %>%
  select(gene_id,`Gene name`) %>%
  filter (! duplicated(gene_id)) %>%
  column_to_rownames("gene_id")

final <- final %>% rownames_to_column("gene_id") %>%
  mutate(gene_id = gsub("\\..*","",.[[1]])) %>%
  left_join(select(genes_match, "Gene name", "Gene stable ID", "Gene description"), by=c("gene_id" = "Gene stable ID")) %>%
  unique() %>% 
  `row.names<-`(., NULL) %>% 
  column_to_rownames("gene_id")

cardio_genes_with_ge <- final %>% rownames_to_column %>%
  filter(abs(log2_Fold_change) > 1 & `p-value(BH_adjusted)` < 0.01 & rowname %in% rownames(cardio_genes)) %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames("rowname")

#########################
########Plotting#########
#########################
ggplot(final, aes(log2_Fold_change, 
                  -log10(`p-value(BH_adjusted)`))) +
  geom_point(data=final[logFCthreshold & pvalueThreshold, ],
             colour="orange", alpha=0.5, size=3) + 
  geom_point(data=final[!logFCthreshold | !pvalueThreshold, ],
             colour="gray", alpha=0.5, size=3) + 
  geom_text_repel(data = final[rownames(cardio_genes_with_ge),], aes(label=`Gene name`),
                 box.padding = 0.4, size=5) +
  theme_light(16) +
  ylab("-log10(|p-value (BH adjusted)|)")

##################################
##########Functional analysis ####
##################################
data=final[logFCthreshold & pvalueThreshold, ]
final_table_ge <- data %>% rownames_to_column("gene_id") %>%
  select(gene_id, `Gene name`, `Gene description`, log2_Fold_change, `p-value(BH_adjusted)`)

source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
fgseaResTidy_hallmarks <- run_fgsea_analysis(final, "~/analysis/genome_utilities/GSEA/h.all.v7.0.symbols.gmt", "Hallmark_genesets", is_deseq = F, is_psichomics = T)
fgseaResTidy_reactome <- run_fgsea_analysis(final, "~/analysis/genome_utilities/GSEA/c2.cp.reactome.v7.0.symbols.gmt", "reactome_genesets", is_deseq = F, is_psichomics = T)
fgseaResTidy_kegg <-  run_fgsea_analysis(final, "~/analysis/genome_utilities/GSEA/c2.cp.kegg.v7.0.symbols.gmt", "kegg_genesets", is_deseq = F, is_psichomics = T)
fgseaResTidy_GO_BP <- run_fgsea_analysis(final, "~/analysis/genome_utilities/GSEA/c5.bp.v7.0.symbols.gmt", "GO_Biological_process", is_deseq = F, is_psichomics = T)
fgseaResTidy_GO_BP <- run_fgsea_analysis(final, "~/analysis/genome_utilities/GSEA/c5..v7.0.symbols.gmt", "GO_Molecular_Function", is_deseq = F, is_psichomics = T)

final_table_ge_up = final_table_ge[final_table_ge$log2_Fold_change > 0,]
final_table_ge_down = final_table_ge[final_table_ge$log2_Fold_change < 0,]
goseq.results <- run_goseq_analysis(final, "hg19", final_table_ge$gene_id, is_deseq = F)
goseq_up_de_genes <- run_goseq_analysis(final, "hg19", final_table_ge_up$gene_id, is_deseq = F)
gpseq_down_de_genes <- run_goseq_analysis(final, "hg19", final_table_ge_down$gene_id, is_deseq = F)
##########################
#####Writing Tables#######
##########################
#%>% filter(gene_id %in% rownames(cardio_genes_with_ge))
write_xlsx(final_table_ge, path  = "GeneExpression_Patients_vs_CTRL.xlsx", col_names=T, format_headers = T)
