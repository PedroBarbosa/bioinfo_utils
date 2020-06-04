library(psichomics)
library(tidyverse)
library(tidylog)
library(writexl)
library(ggplot2)
library(ggrepel)
library(limma)
library(plyr)
library(tidylog)
setwd("/Users/pbarbosa/MEOCloud/analysis/yvan_rna_seq/psichomics/star_counts")
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
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
variance <- apply(psi, 1, var, na.rm = T)
psi <- cbind(psi, variance)
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
event <- "A5SS_14_+_20345087_20345126_20345394_PARP2"
info  <- queryEnsemblByEvent(event, species="human", assembly="hg38")
plotTranscripts(info, event = event)
plotDistribution(psi["A5SS_14_+_20345087_20345126_20345394_PARP2",], psi = T, rug = T, groups = groups)

#####################
###Get coordinates###
#####################
get_coord <- function(x){
  string <- strsplit(x, "_")[[1]]
  event_type <- string[1]
  chrom <- paste0("chr",string[2])
  numbers <- as.integer(string[-c(1,2,3,length(string))])
  coord <- paste0(chrom,":",min(numbers), "-", max(numbers))
  return(c(coord, event_type))
}
strand_map <- list("+" = "plus", "-" = "minus")

#####################
#######Get bed#######
#####################
get_bed <- function(x, comparison){
  coord <- strsplit(x['spanning_coordinate'], ":|-")
  event_type <- ifelse(x["event_type"] == "SE", "ES", x["event_type"])
  return(as_tibble(t(c(coord[[1]][1], coord[[1]][2], coord[[1]][3], paste(x["symbol"], event_type, sep="_"), comparison, x['Strand']))))
}

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
                            & abs(diffSplicing_dcm_ctrl$`∆ Median`) > 0.2)  

useful_info <- lapply(rownames(deltaPSI_dcm_ctrl), get_coord)
deltaPSI_dcm_ctrl$spanning_coordinate <- unlist(map(useful_info, 1))
deltaPSI_dcm_ctrl$event_type <- unlist(map(useful_info, 2))

final_dcm_ctrl <- deltaPSI_dcm_ctrl %>% rownames_to_column("event_id") 
new_symbol <- plyr::mapvalues(final_dcm_ctrl$Gene,from=c("LEPREL2", "FLJ31306", "KIAA0226", "H2AFY", "ZZZ3", "MKL2", "FAM115A", "PPP2R4", "GRAMD3", "GAS7","ATP5J"),
                              to=c("P3H3", "PSMA3-AS1", "RUBCN", "MACROH2A1","ZZZ3","MRTFB","TCAF1","PTPA", "GRAMD2B", "GAS7","ATP5PF"))
final_dcm_ctrl$symbol <- new_symbol
final_dcm_ctrl$is_DCM <- TRUE
final_dcm_ctrl$is_ICM <- FALSE
final_dcm_ctrl <- final_dcm_ctrl %>%left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>% distinct() %>% filter(!is.na(gene_id)) %>%
  dplyr::select(.,c(1, symbol, gene_id, 4, 9, 10, 12, 13, 20, 21, 23, spanning_coordinate, event_type, is_DCM)) %>% 
  write_excel_csv(., "DCM_vs_CTRL.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")

bed <- apply(final_dcm_ctrl, 1, get_bed, "ctrl_dcm")
bed_df_ctrl_dcm <- ldply (bed, data.frame)
bed_df_ctrl_dcm %>%
  dplyr::distinct()%>%
  dplyr::arrange() %>%
  write_tsv("DCM_vs_Ctrl.bed" , col_names = F, na = "NA", append = FALSE, quote_escape = "double")

final_dcm_ctrl$Strand <- sapply(final_dcm_ctrl$Strand, function(x) ifelse(x =="+", "plus","minus"))
final_dcm_ctrl %>% dplyr::select(.,c(spanning_coordinate, symbol, event_type, Strand)) %>%
  write_excel_csv(., "psichomics_to_ggsashimi_DCM_vs_CTRL.csv" , col_names = F, na = "NA", append = FALSE, delim = "\t", quote_escape = "double")
write_xlsx(final_dcm_ctrl, path  = "DCM_vs_CTRL.xlsx", col_names=T, format_headers = T)


#######################
#######ICM vs Ctrl ####
#######################
icm_ctrl <- groups[c("ICM", "Ctrl")]
diffSplicing_icm_ctrl <- diffAnalyses(psi, groups = icm_ctrl, analyses = "wilcoxRankSum")
deltaPSI_icm_ctrl <- subset(diffSplicing_icm_ctrl, diffSplicing_icm_ctrl$`Wilcoxon p-value (BH adjusted)` < 0.1 & 
                              abs(diffSplicing_icm_ctrl$`∆ Median`) > 0.2)   

useful_info <- lapply(rownames(deltaPSI_icm_ctrl), get_coord)
deltaPSI_icm_ctrl$spanning_coordinate <- unlist(map(useful_info, 1))
deltaPSI_icm_ctrl$event_type <- unlist(map(useful_info, 2))

final_icm_ctrl <- deltaPSI_icm_ctrl %>% rownames_to_column("event_id")
new_symbol <- plyr::mapvalues(final_icm_ctrl$Gene,from=c("LEPREL2", "FLJ31306", "KIAA0226", "H2AFY", "ZZZ3", "MKL2", "FAM115A", "PPP2R4", "GRAMD3", "GAS7","ATP5J"),
                              to=c("P3H3", "PSMA3-AS1", "RUBCN", "MACROH2A1","ZZZ3","MRTFB","TCAF1","PTPA", "GRAMD2B", "GAS7","ATP5PF"))
final_icm_ctrl$symbol <- new_symbol
final_icm_ctrl$is_ICM <- TRUE
final_icm_ctrl$is_DCM <- FALSE
final_icm_ctrl <- final_icm_ctrl %>%left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>% distinct() %>% filter(!is.na(gene_id)) %>%
  dplyr::select(.,c(1, symbol,gene_id, 4, 9, 10, 12, 13, 20, 21, 23, spanning_coordinate, event_type, is_ICM)) %>% 
  write_excel_csv(., "ICM_vs_CTRL.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")

bed <- apply(final_icm_ctrl, 1, get_bed, "ctrl_icm")
bed_df_ctrl_icm <- ldply (bed, data.frame)
final_icm_ctrl$Strand <- sapply(final_icm_ctrl$Strand, function(x) ifelse(x =="+", "plus","minus"))
final_icm_ctrl %>% dplyr::select(.,c(spanning_coordinate, Gene, event_type, Strand)) %>%
  write_excel_csv(., "psichomics_to_ggsashimi_ICM_vs_CTRL.csv" , col_names = F, na = "NA", append = FALSE, delim = "\t", quote_escape = "double")
write_xlsx(final_icm_ctrl, path  = "ICM_vs_CTRL.xlsx", col_names=T, format_headers = T)

#######################
#######DCM vs ICM ####
#######################
#MISSING LAST UPDATE
dcm_icm <- groups[c("DCM", "ICM")]
diffSplicing_dcm_icm <- diffAnalyses(psi, groups = dcm_icm, analyses = "wilcoxRankSum")
deltaPSI_dcm_icm <- subset(diffSplicing_dcm_icm, diffSplicing_dcm_icm$`Wilcoxon p-value` < 0.1 & 
                             abs(diffSplicing_dcm_icm$`∆ Median`) > 0.2)   

useful_info <- lapply(rownames(deltaPSI_dcm_icm), get_coord)
deltaPSI_dcm_icm$spanning_coordinate <- unlist(map(useful_info, 1))
deltaPSI_dcm_icm$event_type <- unlist(map(useful_info, 2))

final_dcm_icm <- deltaPSI_dcm_icm %>% rownames_to_column("event_id") %>%
  select(.,c(1, 5, 4, 9, 10, 12, 13, 20, 21, 23, spanning_coordinate, event_type)) %>% 
  arrange_at(ncol(.), desc) %>%
  write_excel_csv(., "DCM_vs_ICM.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")

final_dcm_icm$Strand <- sapply(final_dcm_icm$Strand, function(x) ifelse(x =="+", "plus","minus"))
final_dcm_icm %>% dplyr::select(.,c(spanning_coordinate, Gene, event_type, Strand)) %>%
  write_excel_csv(., "psichomics_to_ggsashimi_DCM_vs_ICM.csv" , col_names = F, na = "NA", append = FALSE, delim = "\t", quote_escape = "double")
write_xlsx(final_dcm_icm, path  = "DCM_vs_ICM.xlsx", col_names=T, format_headers = T)


#######################
#######CONCAT##########
#######################
final_df <- bind_rows(final_dcm_ctrl, final_icm_ctrl) %>% 
  dplyr::select(.,c(1, 2, 3, 4, 13, 7, 8, 12, 14, 15, 9, 10, 17,11)) %>%
  arrange(event_id) %>%
  write_excel_csv(., "all_concat_vs_Ctrl.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")
write_xlsx(final_df, path  = "all_concat_vs_Ctrl.xlsx", col_names=T, format_headers = T)

final_df$Strand <- sapply(final_df$Strand, function(x) ifelse(x =="+", "plus","minus")) 
final_df %>% dplyr::select(.,c(spanning_coordinate, symbol, event_type, Strand)) %>%
  distinct() %>%
  write_excel_csv(., "psichomics_to_ggsashimi.csv" , col_names = F, na = "NA", append = FALSE, delim = "\t", quote_escape = "double")

final_bed <- bed_df_ctrl_dcm %>% full_join(bed_df_ctrl_icm, by=c(colnames(bed_df_ctrl_dcm)[[1]],colnames(bed_df_ctrl_dcm)[[2]],colnames(bed_df_ctrl_dcm)[[3]]))
final_bed[] <- t(apply(final_bed, 1, function(x) `length<-`(na.omit(x), length(x)))) #shift NAs to the left
final_bed <- unite(final_bed, groups, c(colnames(final_bed)[[5]], colnames(final_bed)[[8]]), sep=";") %>%
  mutate_all(~gsub(";NA","",.)) %>%
  dplyr::select(c(1:6)) %>%
  distinct() %>%
  arrange(., .[[1]], .[[2]]) %>%
  write_tsv("all_concat_vs_Ctrl.bed" , col_names = F, na = "NA", append = FALSE, quote_escape = "double")

#######################
####GENE EXPRESSION####
#######################
geneExprFiltered <- geneExpr[rowSums(geneExpr)>100,]
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
####PATIENTS VS CONTROLS######
to_match <- c("DCM", "ICM")
groups_sick_vs_ctrl <- list("Ctrl"=grep("Ctrl", colnames(geneExpr), value=TRUE), 
               "Patient"= grep("Ctrl", colnames(geneExpr), invert = T, value=TRUE))

pgenes <- rownames(geneExprNorm)
plotDistribution(geneExprNorm["ENSG00000163359", ], groups, psi=FALSE)

groups_dcm_vs_icm <- list("DCM"=grep("DCM", colnames(geneExpr), value = TRUE),
                      "ICM"=grep("ICM", colnames(geneExpr), value = TRUE))

ge_dcm_icm <- geneExprNorm[, unlist(groups_dcm_vs_icm), drop=FALSE]
ge_sick_ctrl <- geneExprNorm[, unlist(groups_sick_vs_ctrl), drop=FALSE]

isFromGroup1 <- colnames(ge_sick_ctrl) %in% groups_sick_vs_ctrl[[1]]
isFromGroup1 <- colnames(ge_dcm_icm) %in% groups_dcm_vs_icm[[1]]

design <- cbind(1, ifelse(isFromGroup1, 0, 1))
fit <- lmFit(ge_sick_ctrl, design = design)
fit <- lmFit(ge_dcm_icm, design = design)

ebayes_fit <- eBayes(fit, trend=T)
pvalueAdj <- "BH"
summary <- topTable(ebayes_fit, number=nrow(fit), coef=2, sort.by = "none", adjust.method = pvalueAdj, confint = T)
names(summary) <- c("log2_Fold_change", "conf_int1", "conf_int2", "moderated_t-statistics", "avg_expression","p-value", 
                    paste0("p-value(",pvalueAdj,"_adjusted)"), "B-statistics")
attr(summary, "groups") <- groups_sick_vs_ctrl
attr(summary, "groups") <- groups_dcm_vs_icm
stats <- diffAnalyses(ge_sick_ctrl, groups_sick_vs_ctrl, "basicStats", geneExpr = T, pvalueAdjust = NULL)
stats <- diffAnalyses(ge_dcm_icm, groups_dcm_vs_icm, "basicStats", geneExpr = T, pvalueAdjust = NULL)
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
