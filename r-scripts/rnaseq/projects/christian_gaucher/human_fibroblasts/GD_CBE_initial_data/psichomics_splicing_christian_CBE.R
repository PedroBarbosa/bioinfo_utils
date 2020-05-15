library(psichomics)
library(tidyverse)
library(tidylog)
library(writexl)
library(ggplot2)
library(ggrepel)
library(limma)
#################################
#########Preparing the data #####
#################################
setwd("/Users/pbarbosa/MEOCloud/analysis/christian/standard_rna_seq_analysis/human/splicing/psichomics")

#fixInNamespace("prepareGeneQuantSTAR", pos="package:psichomics")
#setnames(table, colnames(table)[[2]], paste0("col", index))
prepareGeneQuant("GD1"="GD1ReadsPerGene.out.tab",
                 "GD2"="GD2ReadsPerGene.out.tab",
                 "GD3"="GD3ReadsPerGene.out.tab",
                 "S_3y_CBE"="S_3y_CBEReadsPerGene.out.tab",
                 "S_10y_CBE"="S_10y_CBEReadsPerGene.out.tab",
                 "S_11y_CBE"="S_11y_CBEReadsPerGene.out.tab",
                 "S_3y_CTRL"="S_3y_CTRLReadsPerGene.out.tab",
                 "S_10y_CTRL"="S_10y_CTRLReadsPerGene.out.tab",
                 "S_11y_CTRL"="S_11y_CTRLReadsPerGene.out.tab",
                 strandedness = "stranded (reverse)")


prepareJunctionQuant("GD1SJ.out.tab",
                     "GD2SJ.out.tab",
                     "GD3SJ.out.tab",
                     "S_3y_CBESJ.out.tab",
                     "S_10y_CBESJ.out.tab",
                     "S_11y_CBESJ.out.tab",
                     "S_3y_CTRLSJ.out.tab",
                     "S_10y_CTRLSJ.out.tab",
                     "S_11y_CTRLSJ.out.tab"
                     )

#####################################
#########Loading the data ###########
#####################################
data <- loadLocalFiles("/Users/pbarbosa/MEOCloud/analysis/christian/standard_rna_seq_analysis/human/splicing/psichomics/")
geneExpr <- data$Data$`Gene expression`
junctionQuant <- data$Data$`Junction quantification`

groups <- list("GD"=colnames(geneExpr)[1:3], 
               "CBE"=colnames(geneExpr)[4:6], 
                "CTRL"=colnames(geneExpr)[7:9])

human <- listSplicingAnnotations()[[3]]
annotation <- loadAnnotation(human)

#####################
####Preprocessing####
#####################
psi <- quantifySplicing(annotation, junctionQuant, minReads=10)
events <- rownames(psi)

psi <- psi[rowSums(is.na(psi)) <= 1,]

variance <- apply(psi, 1, var, na.rm = T)
psi <- cbind(psi, variance)

psi <- psi[ psi$variance > quantile(psi$variance , 0.90, na.rm = T) , ]
psi <- psi[, !(colnames(psi) %in% c("variance"))]
psi_pca <- performPCA(t(psi))

#manual pca
percentage <- round(psi_pca$sdev / sum(psi_pca$sdev) * 100, 2)
psi_pca_manual <- as.data.frame(psi_pca$x)
group <- c(rep("GD", 3), rep("CBE", 3), rep("CTRL", 3))
psi_pca_manual$group <- group
percentage <- paste( colnames(psi_pca_manual), "(", paste0( as.character(percentage), "%", ")", sep="") )

p <- ggplot(psi_pca_manual,aes(x=PC1, y=PC2, color=group))
p<-p+geom_point(size=3) + xlab(percentage[1]) + ylab(percentage[2])
p


######################
######ALL GROUPS######
######################
diffSplicing_all <- diffAnalyses(psi, groups = groups, analyses = "kruskal")
deltaPSIthreshold_all <- diffSplicing_all$`Kruskal p-value (BH adjusted)` < 0.05
delta_psi_all <- diffSplicing_all[which(deltaPSIthreshold_all), ]
dim(delta_psi_all)


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
  coord = strsplit(x['spanning_coordinate'], ":|-")
  return(as_tibble(t(c(coord[[1]][1], coord[[1]][2], coord[[1]][3], paste(x["Gene"], x["Event type"], sep="_"), comparison, x['Strand']))))
}


######################
#####CBE vs Ctrl #####
######################
CBE_ctrl <- groups[c("CBE", "CTRL")]
diffSplicing_CBE_ctrl <- diffAnalyses(psi, groups = CBE_ctrl, pvalueAdjust= "BH", analyses = "wilcoxRankSum")

deltaPSI_CBE_ctrl <- subset(diffSplicing_CBE_ctrl, diffSplicing_CBE_ctrl$`Wilcoxon p-value` < 0.1
                           & abs(diffSplicing_CBE_ctrl$`∆ Median`) > 0.2)  

useful_info <- lapply(rownames(deltaPSI_CBE_ctrl), get_coord)
deltaPSI_CBE_ctrl$spanning_coordinate <- unlist(map(useful_info, 1))
deltaPSI_CBE_ctrl$event_type <- unlist(map(useful_info, 2))
 
final_cbe_ctrl <- deltaPSI_CBE_ctrl %>% rownames_to_column("event_id") %>%
  select(.,c(5, 4, 25, 9, 10, 12, 13, 20, 21, 23, 24)) %>% 
  arrange_at(ncol(.), desc) %>%
  write_excel_csv(., "CBE_vs_CTRL_psichomics.csv" , na = "NA", append = FALSE,
                  delim = "\t", quote_escape = "double")
 write_xlsx(final_cbe_ctrl, path  = "CBE_vs_CTRL_psichomics.xlsx", col_names=T, format_headers = T)
plotDistribution(psi["A3SS_16_-_31094557_31093321_31093482_PRSS53/VKORC1", ], groups) 


bed <- apply(final_cbe_ctrl, 1, get_bed, "CBE_vs_mock")
bed_df_cbe_mock <- ldply (bed, data.frame)
bed_df_cbe_mock %>%
  dplyr::distinct()%>%
  dplyr::arrange() %>%
  write_tsv("CBE_vs_mock.bed" , col_names = F, na = "NA", append = FALSE, quote_escape = "double")

final_cbe_ctrl$Strand <- sapply(final_cbe_ctrl$Strand, function(x) ifelse(x =="+", "plus", "minus"))
final_cbe_ctrl %>% dplyr::select(.,c(spanning_coordinate, Gene, event_type, Strand)) %>%
  write_excel_csv(., "psichomics_to_ggsashimi_CBE_vs_mock.csv" , col_names = F, na = "NA", append = FALSE, delim = "\t", quote_escape = "double")


#######################
####GENE EXPRESSION####
#######################
geneExprFiltered <- geneExpr[rowSums(geneExpr)>50,]
geneExprNorm <- normaliseGeneExpression(geneExprFiltered, log2transform=TRUE)
pcaGE_all <- performPCA(t(geneExprNorm))

#Set gene ID as rownames
rownames(final_df) <- df[,1]
final_df[,1] <- NULL


###################
####manual pca#####
###################
ge_pca <- prcomp(t(geneExprNorm), center = T,scale. = T)
percentage <- round(ge_pca$sdev / sum(ge_pca$sdev) * 100, 2)
ge_pca <- as.data.frame(ge_pca$x)
ge_pca$group <- group
percentage <- paste( colnames(ge_pca), "(", paste0( as.character(percentage), "%", ")", sep="") )

p <- ggplot(ge_pca,aes(x=PC1, y=PC2, color=group))
p<-p+geom_point(size=3) + xlab(percentage[1]) + ylab(percentage[2])
p

#############
####groups###
#############
####PATIENTS VS CONTROLS######
to_match <- c("GD", "CTRL")

pgenes <- rownames(geneExprNorm)
groups_cbe_ctrl <- list("GD"=grep("GD", colnames(geneExpr), value = TRUE),
                      "CTRL"=grep("CTRL", colnames(geneExpr), value = TRUE))

ge_cbe_ctrl <- geneExprNorm[, unlist(groups_cbe_ctrl), drop=FALSE]

isFromGroup1 <- colnames(ge_cbe_ctrl) %in% groups_cbe_ctrl[[1]]
design <- cbind(1, ifelse(isFromGroup1, 0, 1))
fit <- lmFit(ge_cbe_ctrl, design = design)
ebayes_fit <- eBayes(fit, trend=T)
pvalueAdj <- "BH"
summary <- topTable(ebayes_fit, number=nrow(fit), coef=2, sort.by = "none", adjust.method = pvalueAdj, confint = T)
names(summary) <- c("log2_Fold_change", "conf_int1", "conf_int2", "moderated_t-statistics", "avg_expression","p-value", 
                    paste0("p-value(",pvalueAdj,"_adjusted)"), "B-statistics")
attr(summary, "groups") <- groups_cbe_ctrl

stats <- diffAnalyses(ge_cbe_ctrl, groups_cbe_ctrl, "basicStats", geneExpr = T, pvalueAdjust = NULL)
final <- cbind(stats, summary)


#Statistical significance
logFCthreshold  <- abs(final$`log2_Fold_change`) > 1
pvalueThreshold <- final$`p-value(BH_adjusted)` < 0.01


#########################
########Plotting#########
#########################
ggplot(final, aes(log2_Fold_change, 
                  -log10(`p-value(BH_adjusted)`))) +
  geom_point(data=final[logFCthreshold & pvalueThreshold, ],
             colour="orange", alpha=0.5, size=3) + 
  geom_point(data=final[!logFCthreshold | !pvalueThreshold, ],
             colour="gray", alpha=0.5, size=3) + 
  theme_light(16) +
  ylab("-log10(|p-value (BH adjusted)|)")

