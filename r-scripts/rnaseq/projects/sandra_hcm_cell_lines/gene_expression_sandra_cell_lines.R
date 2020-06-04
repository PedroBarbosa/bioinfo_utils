library(edgeR)
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(biomaRt)
library(RColorBrewer)
library(ggrepel)
library(tximport)
library(gtools)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
################################################
########### CARDIOMYOCYTES CELL LINES ##########
################################################
#####################
#####From Salmon#####
#####################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/sandra_HCM_rnaseq/salmon/")
files <- list.files(".", pattern = "*sf")
names(files) <- sapply(strsplit(files,"\\."), function(x) x[1])


tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt")
tx_ids <- readr::read_tsv(file = "DF6.sf") %>% .$Name 
tx_ids <- sapply(strsplit(tx_ids,"\\."), function(x) x[1])
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID`, `Gene stable ID`) , by = c("value" = "Transcript stable ID")) 
txi_salmon <- tximport(files, type = "salmon", txIn = T, txOut = F, tx2gene = tx2gene, ignoreTxVersion = T)

group <- as.factor(rep(c("ctrl", "xutl"), times = c(3,1)))
coldata <- data.frame(group,row.names=colnames(txi_salmon$counts))

dds <- DESeqDataSetFromTximport(txi = txi_salmon,
                         colData = coldata,
                         design = ~ group)
keep <- rowSums(counts(dds)) >= 40 
dds <- dds[keep,]

#########################
#####DATA EXPLORATION####
#########################
group_combination <- c("group","xutl", "ctrl")
run_analysis(dds, group_combination, 1, 0.05, "hg38", explore_data = T)

rld <- rlog(dds)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))
percentVar <- pc$sdev^2 / sum( pc$sdev^2 )

pca_data <- DESeq2::plotPCA(rld, intgroup = "group", returnData=T)
pca_matrix <- as.data.frame(pc$x)
pca_matrix$group <- pca_data$group
ggplot(pca_matrix,aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=5) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  theme(legend.title=element_blank(), text = element_text(size = 20))


################################
####### CORRELATION PLOTS ######
################################
rld <- as.data.frame(assay(rld))
combinations <- as_tibble(combinations(ncol(rld), 2, colnames(rld)))

plotlist=list()
i=1
for(i in 1:nrow(combinations)) {
  pair1 <- as.character(combinations[i, 1])
  pair2 <- as.character(combinations[i, 2])
  
  #For RNA-seq counts, however, the expected variance grows with the mean. 
  #For example, if one performs PCA directly on a matrix of counts
  #or normalized counts (e.g. correcting for differences in sequencing depth), 
  #the resulting plot typically depends mostly on the genes with highest counts because 
  #they show the largest absolute differences between samples.
  #A simple and often used strategy to avoid this is to take the logarithm of the normalized count values
  #plus a pseudocount of 1; however, depending on the choice of pseudocount, 
  #now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot,
  #because taking the logarithm of small counts actually inflates their variance. 
  #p<-ggplot(norm_counts_salmon, aes(x=log2(norm_counts_salmon[,pair1]), y=log2(norm_counts_salmon[, pair2]))) +
  p<-ggplot(rld, aes(x=rld[,pair1], y=rld[, pair2])) +
    geom_point(color='brown', alpha=0.7) +
    labs(x=pair1, y=pair2) +
    stat_cor(method="pearson", digits=3) + 
    geom_smooth(method = "lm", colour="black", size=0.5) +
    theme(text = element_text(size =15))
  
  plotlist[[i]] = ggplotGrob(p)
  i= i + 1 
} 
grid.arrange(grobs=plotlist,nrow=2,ncol=3,
             #bottom="Pairwise comparisons of gene expression (rlog transformation of count data)", 
             gp=gpar(fontface="bold", col="YELLOW", fontsize=20))




#################################
######### DE analysis ###########
#################################
dds_norm <- DESeq(dds)
group_combination <- c("group","xutl", "ctrl")
log2cutoff = 1
padjcutoff = 0.05

list_de_salmon <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, use_contrast = F, FC_shrinkage_method = "apeglm", 
                               genome = "hg38", annotate_locally = T, explore_data = F)
de_genes <- list_de_salmon[[2]]
all_annot_salmon <- annotate_results(list_de_salmon[[1]], T, "hg38") 
all_annot_salmon[all_annot_salmon$symbol == "MYH7",]

#########################################
######## PLOT individual expression #####
#########################################
target_genes <- c("MYH7", "MYBPC3", "TNNT2", "TNNI3", "ACTC", "TPM1", "MYL3", "MYL2", "TTN", "MYH6")
target_genes <- c("TBC1D3")
rownames(dds_norm) <- gsub("\\..*", "", rownames(dds_norm))
norm_counts <- as.data.frame(counts(dds_norm, normalize=T)) %>% rownames_to_column(var="gene_id") %>%
  left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>% distinct()

counts_target_genes <- norm_counts %>% filter(symbol %in% target_genes)

long_table <-  pivot_longer(counts_target_genes, 
                            cols = colnames(counts_target_genes)[! colnames(counts_target_genes) %in% c("gene_id", "symbol")], 
                            values_to = "norm_expression", names_to = "sample")
long_table$group <- ifelse(long_table$sample == "Xutl", "Xutl", "Ctrl")

out <- long_table %>% group_by(gene_id) %>%
  do(plots=ggplot(data =., aes(x=sample, y=norm_expression, fill = group)) +
       geom_bar(stat="identity",alpha=0.7, color = "black") +
       ggtitle(.$symbol) +
       labs(x = "Cell line", y = "Norm expression (Median of ratios)") +
       theme(text = element_text(size=15), aspect.ratio = 1.2))

out$plots[[1]]
###################################
####### Expression in chrX ########
###################################
tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt")
tpm_gene_expression = as.data.frame(txi_salmon$abundance)

tpm_gene_expression_chrX <- rownames_to_column(tpm_gene_expression, var="gene_id") %>%
        left_join(dplyr::select(tx2gene, `Gene stable ID`, `Chromosome/scaffold name`, `Gene name`), by = c("gene_id" = "Gene stable ID")) %>%
        rename(chrom = `Chromosome/scaffold name`, gene_name = `Gene name`) %>% distinct() %>% 
        filter(chrom == "X") %>% 
        dplyr::select(-c(Xutl, chrom))

TCLAB <- tpm_gene_expression_chrX %>% filter(TCLAB > 50) %>% dplyr::select(gene_id, gene_name, TCLAB) %>% arrange(-TCLAB) 
#%>%
#  write_xlsx(., path = "~/Desktop/TCLAB.xlsx", col_names=T, format_headers = T)
DF6 <- tpm_gene_expression_chrX %>% filter(DF6 > 50) %>% dplyr::select(gene_id, gene_name, DF6) %>% arrange(-DF6) 
#%>%
#  write_xlsx(., path = "~/Desktop/DF6.xlsx", col_names=T, format_headers = T)
GIBCO <- tpm_gene_expression_chrX %>% filter(GIBCO > 50) %>% dplyr::select(gene_id, gene_name, GIBCO) %>% arrange(-GIBCO) 
#%>%
#  write_xlsx(., path = "~/Desktop/GIBCO.xlsx", col_names=T, format_headers = T)

isec <- Reduce(intersect, list(TCLAB$gene_name, DF6$gene_name, GIBCO$gene_name)) 
write(isec, file = "~/Desktop/isec.txt")

