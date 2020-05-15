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
#######################################
###########CBE vs CTRL PAPER ##########
#######################################
#####################
#FROM FEATURE COUNTS#
#####################
feature_counts_data = read.table("~/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/christian/human_fibroblasts/cbe_ctrl_april_2019/star/uniq/feature_counts/featureCounts.txt",
                                 header=T, sep="\t")

#n_occur <- data.frame(table(features_counts_data$`gene_name`))
#print(n_occur[n_occur$Freq > 1,])

#Remove gene_name col
df_numbers = feature_counts_data[ , -which(names(feature_counts_data) %in% c("gene_name", "Chr", "Start", "End", "Strand", "Length"))]
colnames(df_numbers) <- gsub("X.home.pedro.barbosa.mcfonseca.shared.christian.human_fibroblasts.cbe_ctrl_april_2019.star.uniq.", "", colnames(df_numbers)) 
colnames(df_numbers) <- gsub("_uniq.bam", "", colnames(df_numbers))

#Set gene id as rowname and remove GD columns
df_numbers <- column_to_rownames(df_numbers, "Geneid")
df_numbers <- df_numbers[, !grepl("GD", colnames(df_numbers))]

treatment <- as.factor(rep(c("CBE", "mock"), times = 3))
cellline <- as.factor(rep(c("F1", "F2", "F3"), each = 2))
coldata <- data.frame(treatment, cellline,row.names=colnames(df_numbers))
dds_fc <- DESeqDataSetFromMatrix(countData = df_numbers,
                                 colData = coldata,
                                 design = ~ cellline + treatment)

keep <- rowSums(counts(dds_fc) >= 10) >= 4
dds_fc <- dds_fc[keep,]
dds_fc$treatment <- relevel(dds_fc$treatment,ref="mock")
group_combination <- c("treatment","CBE", "mock")

#####################
#####From Salmon#####
#####################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/christian/human_fibroblasts/cbe_ctrl_april_2019/salmon")
files <- list.files(".", pattern = "*sf")

treatment <- as.factor(unlist(lapply(strsplit(files, split = "_"), `[`, 2)))
cellline <-  as.factor(unlist(lapply(strsplit(files, split = "_"), `[`, 1)))
samples <-  as.factor(unlist(lapply(strsplit(files, split = "_quant"), `[`, 1)))
names(files) <- names(files) <- samples
coldata <- data.frame(treatment, cellline, row.names = samples)

tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt")
tx_ids <- readr::read_tsv(file = "F1_CBE_quant.sf") %>% .$Name 
tx_ids <- sapply(strsplit(tx_ids,"\\."), function(x) x[1])
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID`, `Gene stable ID`) , by = c("value" = "Transcript stable ID")) 
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene, ignoreTxVersion = T)

#Multi Factor design to include sample information as a term
#This will account for differences between the samples while estimating the effect due to the treeatment
dds_salmon <- DESeqDataSetFromTximport(txi = txi_salmon,
                                       colData = coldata,
                                       design = ~ cellline + treatment)

keep <- rowSums(counts(dds_salmon) >= 10) >= 4
#keep <- rowSums(counts(dds_salmon)) >= 10
dds_salmon <- dds_salmon[keep,]
dds_salmon$treatment <- relevel(dds_salmon$treatment,ref="mock")
group_combination <- c("treatment","CBE", "mock")


#########################
#####DATA EXPLORATION####
#########################
run_analysis(dds_fc, group_combination, 1, 0.05, "hg38", explore_data = T)
#Dendogram
dists <- dist(t(assay(rld_fc)))
plot(hclust(dists))

########################
### 2D PCA from DESeq ##
########################
#dds_norm_salmon <- DESeq(dds_salmon)
#norm_counts_salmon <- 1 + as.data.frame(counts(dds_norm_salmon, normalized=T))
#dds_norm_fc <- DESeq(dds_fc)
#norm_counts_fc <- 1 + as.data.frame(counts(dds_norm_fc, normalized=F))

rld_salmon <- rlog(dds_salmon)
rld_fc <- rlog(dds_fc)
DESeq2::plotPCA(rld_fc, ntop=500, intgroup = "treatment")
pca_data <- DESeq2::plotPCA(rld_fc, intgroup = "treatment", returnData=T)

ggplot(pca_data,aes(x=PC1, y=PC2, color=name)) +
  geom_point(size=5) + 
  scale_color_manual(values=c("skyblue3", "skyblue4", "darkseagreen3", "darkseagreen4", "brown3", "brown4"))

########################
##### Manual 2DPCA #####
########################
rv <- rowVars(assay(rld_salmon))
rv <- rowVars(assay(rld_fc))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(rld_fc)[select,]))
percentVar <- pc$sdev^2 / sum( pc$sdev^2 )

pca_matrix <- as.data.frame(pc$x)
pca_matrix$treatment <- pca_data$treatment
colors <- c("blue", "red") 
colors <- colors[as.numeric(pca_data$treatment)] 
ggplot(pca_matrix,aes(x=PC1, y=PC2, color=rownames(pca_matrix))) +
  geom_point(size=5) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("skyblue3", "skyblue4", "darkseagreen3", "darkseagreen4", "brown3", "brown4")) +
  theme(legend.title=element_blank(), text = element_text(size = 20))

######################
####### 3D PCA #######
######################
library("scatterplot3d")
scatter3D <-scatterplot3d(pc$x[,1:3], xlab="PC1",ylab="PC2", zlab="PC3", grid=F, color = colors,  pch = 16)
legend("right",legend = levels(pca_matrix$groups), col = colors, pch = 16,bty="n")
text(1,  1, labels = rownames(pca_matrix)[1],cex= 0.7, col = "black",pos=2.5)

###############
##Interactive##
###############
library(rgl)
library(knitr)
options(rgl.printRglwidget = TRUE)
options(rgl.useNULL = TRUE)
plot3d(pc$x[,1:3], col=colors, type='s', size = 3, surface=FALSE, ellipsoid = TRUE)
text3d(pc$x[,1] + 2, pc$x[,2] + 2, pc$x[,3] + 2 , texts =rownames(pc$x), adj = 0.1, 
       color="black", family="serif", font=7, cex=1)
setwd("~/Desktop/")
writeASY(outtype = "pdf")
#legend3d("topright", legend = levels(pca_matrix$groups), pch = 16, col = colors, cex=1, inset=c(0.02))
a


###########################
####Correlation plot ######
###########################
##Multiple plot
rld_fc <- as.data.frame(assay(rld_fc))
combinations_mock <- as_tibble(combinations(ncol(rld_fc)/2, 2, colnames(rld_fc)[grepl("mock", colnames(rld_fc))]))
combination_CBE <- as_tibble(combinations(ncol(rld_fc)/2, 2, colnames(rld_fc)[grepl("CBE", colnames(rld_fc))]))
combinations <- rbind(combinations_mock, combination_CBE)

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
  p<-ggplot(rld_fc, aes(x=rld_fc[,pair1], y=rld_fc[, pair2])) +
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


#######################
######### DE ##########
#######################
log2cutoff <- 1
padjcutoff <- 0.05
setwd("~/Desktop/")
#################
#### salmon #####
#################
list_de_salmon <- run_analysis(dds_salmon, group_combination, log2cutoff, padjcutoff, use_contrast = F, FC_shrinkage_method = "apeglm", 
                               genome = "hg38", annotate_locally = T, explore_data = F)

de_genes_salmon_apeglm_shrinkage <- rownames(list_de_salmon[[2]])
all_annot_salmon <- annotate_results(list_de_salmon[[1]], T, "hg38") 

as.data.frame(list_de_salmon[[1]])%>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000117399'))
as.data.frame(counts(dds_salmon), normalized=T) %>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000117399'))


###################
#### FeatureC #####
###################
list_de_fc <- run_analysis(dds_fc, group_combination, log2cutoff, padjcutoff, use_contrast = F, FC_shrinkage_method = "apeglm", 
                           genome = "hg38", annotate_locally = T, explore_data = F)
de_genes_fc_apeglm_shrinkage <- rownames(list_de_fc[[2]])
all_annot_fc <- annotate_results(list_de_fc[[1]], T, "hg38") 

as.data.frame(list_de_fc[[1]])%>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000162383'))
as.data.frame(counts(dds_fc), normalized=T) %>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000162383'))


###OLD ANALYSIS GENE##
library(readxl)
old_analysis <- read_excel("/Users/pbarbosa/MEOCloud/analysis/christian/standard_rna_seq_analysis/human/CBE_vs_CTRL_DESeq2.xlsx") 
genes_old_analysis <- old_analysis %>% pull(gene_id)

##########################################
######## VENN PLOT between DE genes ######
##########################################
library(VennDiagram)
grid.newpage()    
plt <- venn.diagram(x = list(genes_old_analysis, de_genes_fc_apeglm_shrinkage),
                    #x = list(genes_old_analysis, de_genes_fc_apeglm_shrinkage, de_genes_salmon_apeglm_shrinkage)
                    na = "remove",
                    category.names = c("Old analysis", "New analysis "),
                    fill = c('lightyellow', 'paleturquoise'),
                    alpha = c(0.5, 0.5),
                    #category.names = c("Old FC" , "FC normal shrinkage" , "Salmon normal shrinkage"),
                    #fill = c('lightyellow', 'paleturquoise', 'lightgreen'),
                    #alpha = c(0.5, 0.5, 0.5),
                    print.mode = c("raw", "percent"),
                    filename = NULL,
                    output = T)
grid::grid.draw(plt)

#########################################
###### UNIQUE GENES TO EACH APPROACH ####
#########################################
unique_genes_old <- Reduce(setdiff, list(genes_old_analysis, de_genes_fc_apeglm_shrinkage, de_genes_salmon_apeglm_shrinkage))
old_analysis %>% filter(gene_id %in% unique_genes_old) %>%
  write_xlsx(path = paste0("~/Desktop/unique_old_fc",".xlsx", sep=""), 
             col_names=T, format_headers = T)

unique_genes_fc <- Reduce(setdiff, list(de_genes_fc_apeglm_shrinkage, genes_old_analysis, de_genes_salmon_apeglm_shrinkage)) 
unique_genes_fc <- list_de_fc[[2]][unique_genes_fc,] %>%
  dplyr::select(.,c(symbol, description,log2FoldChange)) %>% rownames_to_column(var="gene_id") %>%
  write_xlsx(., path = paste0("~/Desktop/fc_apeglm_shrinkage",".xlsx", sep=""), 
             col_names=T, format_headers = T)

unique_genes_salmon <- Reduce(setdiff, list(de_genes_salmon_apeglm_shrinkage, genes_old_analysis, de_genes_fc_apeglm_shrinkage))
unique_genes_salmon <- list_de_salmon[[2]][unique_genes_salmon,] %>%
  dplyr::select(.,c(symbol, description,log2FoldChange)) %>% rownames_to_column(var="gene_id") %>%
  write_xlsx(., path = paste0("~/Desktop/salmon_apeglm_shrinkage",".xlsx", sep=""), 
             col_names=T, format_headers = T)


#####################################
######## Functional enrichment ######
#####################################
#GPROFILER
mock_cbe_gprofiler <- run_enrichment_gprofiler(de_genes_fc_apeglm_shrinkage,custom_genes_background = rownames(all_annot_fc),
                                               retrieve_only_sign_results = T, 
                                               exclude_iea = F, 
                                               retrieve_short_link = F,
                                               measure_under = T, 
                                               domain_scope = "annotated",
                                               sources= NULL)
dim(mock_cbe_gprofiler)
plot_enrichment_results(mock_cbe_gprofiler, top_n = 20, short_term_size = F, add_info_to_labels = F, size_of_short_term = 500)

#GOseq
goseq.results <- run_goseq_analysis(de_genes_fc_apeglm_shrinkage, rownames(all_annot_fc), "hg19")
plot_enrichment_results(dplyr::filter(out_df_goseq, significant == T), 
                        rank_by = "log10padj_over", top_n = 20, short_term_size = F, size_of_short_term = 200)

#FGSEA
out <- run_fgsea_analysis(all_annot_fc, is_ranked_already = F, top_n = 20, npermutations = 2000)
plot_enrichment_results(dplyr::filter(out, significant == T), 
                        rank_by = "NES", label = "Normalized Enrichment Score", top_n = 20, 
                        font_size = 15,
                        reverse_order = T,
                        add_info_to_labels =  F,
                        short_term_size = F, 
                        size_of_short_term = 50)

plot_enrichment_results(dplyr::filter(out, significant == T & source == "KEGG"), 
                        top_n = 20, short_term_size = F, size_of_short_term = 500)



######################################################################
############## Heatmap of DE genes and overall PCA ###################
######################################################################
rld_fc <- rlog(dds_fc)
rownames(rld_fc) <- gsub("\\..*", "", rownames(rld_fc))
rld_de_genes <- assay(rld_fc[de_genes_fc_apeglm_shrinkage,])
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
plot_expression_heatmaps(rld_de_genes)

dev.off()

#####################################
#######SPLICING INTERSECTION ########
#####################################
rmats_positive_genes <- read_tsv("rmats/rmats_positive_list_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

majiq_positive_genes <- read_tsv("majiq/mock_vs_CBE_majiq_positive_list_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

vast_positive_genes <- read_tsv("vast-tools/2_positive_geneIDs_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

all_positive_genes <- Reduce(dplyr::union, list(vast_positive_genes, majiq_positive_genes, rmats_positive_genes))
write(all_positive_genes, file = "1_splicing_positive_genes_to_GO.txt")
all_positive_genes <- read_tsv("1_splicing_positive_from_pipeline.txt", col_names = F) %>% pull()

dplyr::intersect(all_positive_genes, de_genes_fc_apeglm_shrinkage)

library(VennDiagram)
grid.newpage()    
plt <- venn.diagram(x = list(all_positive_genes, de_genes_fc_apeglm_shrinkage),
                    #x = list(genes_old_analysis, de_genes_fc_apeglm_shrinkage, de_genes_salmon_apeglm_shrinkage)
                    na = "remove",
                    category.names = c("Genes with splicing events", "DE genes"),
                    fill = c('lightyellow', 'paleturquoise'),
                    alpha = c(0.5, 0.5),
                    #category.names = c("Old FC" , "FC normal shrinkage" , "Salmon normal shrinkage"),
                    #fill = c('lightyellow', 'paleturquoise', 'lightgreen'),
                    #alpha = c(0.5, 0.5, 0.5),
                    print.mode = c("raw"),
                    filename = NULL,
                    output = T)
grid::grid.draw(plt)


#############################################
#### DE with Splicing (RI, A3SS or A5SS)#####
#############################################
de_and_splicing <- read_tsv("3_DE_genes_with_splicing.txt", col_names = "geneid") %>% 
  left_join(dplyr::select(rownames_to_column(list_de_fc[[2]], var="geneid"), c(symbol,geneid))) %>% dplyr::select(symbol)
            
target_genes <- read_tsv("5_DE_splicing_concat.tsv", col_names=c("symbol","event_tyoe","dpsi")) %>% 
  left_join(rownames_to_column(list_de_fc[[2]])) %>%
  dplyr::select(symbol, rowname, log2FoldChange) %>% distinct() %>% column_to_rownames()

dds_fc_norm <- DESeq(dds_fc)
rownames(dds_fc_norm) <- gsub("\\..*", "", rownames(dds_fc_norm))
norm_counts<- as.data.frame(counts(dds_fc_norm,normalize=T)[rownames(target_genes),])
norm_counts$mock <- rowMeans(subset(norm_counts, select = grepl("mock", colnames(norm_counts))), na.rm = FALSE)
norm_counts$CBE <- rowMeans(subset(norm_counts, select = grepl("CBE", colnames(norm_counts))), na.rm = FALSE)

long_table <- rownames_to_column(norm_counts,var="geneid") %>% left_join(rownames_to_column(target_genes,var="geneid")) %>%
  pivot_longer(cols = starts_with("F"), values_to = "avg_norm_expression", names_to = "sample") %>%
  dplyr::select(c(symbol, sample, avg_norm_expression)) %>% separate(col = sample, into= c("ind", "group"))
long_table$group <-factor(long_table$group, levels=c("mock", "CBE"))

ggplot(long_table,  aes(x=symbol, y=avg_norm_expression, fill=group)) +
  geom_boxplot() +
  xlab("Gene name") +
  theme(text = element_text(size=20))
  