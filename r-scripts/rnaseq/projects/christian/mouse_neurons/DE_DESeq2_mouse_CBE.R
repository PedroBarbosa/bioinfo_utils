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


source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
#####################
##### From FC #######
#####################
# feature_counts_data = read.table("~/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/christian/mouse_neurons_may_2020/star/uniq/feature_counts/featureCounts.txt",
#                                  header=T, sep="\t")
# 
# #Remove gene_name col
# df_numbers = feature_counts_data[ , -which(names(feature_counts_data) %in% c("gene_name", "Chr", "Start", "End", "Strand", "Length"))]
# colnames(df_numbers) <- gsub("X.home.pedro.barbosa.mcfonseca.shared.christian.mouse_neurons_may_2020.star.uniq.", "", colnames(df_numbers)) 
# colnames(df_numbers) <- gsub("_uniq.bam", "", colnames(df_numbers))
# 
# #Set gene id as rowname and remove GD columns
# df_numbers <- column_to_rownames(df_numbers, "Geneid")
# df_numbers <- df_numbers[, !grepl("NSC", colnames(df_numbers))]
# colnames(df_numbers) <- c("B61_CBE", "B62_CBE", "B63_CBE", "B61_mock", "B62_mock", "B63_mock",
#                           "BC1_CBE", "BC2_CBE", "BC3_CBE", "BC1_mock", "BC2_mock", "BC3_mock")
# 
# treatment <- rep(as.factor(rep(c("CBE", "mock"), times = 3)), times=2)
# cellline <- as.factor(rep(c("B6", "BC"), each = 6))
# replicate <- as.factor(rep(c("1", "2", "3"), times=4))
# coldata <- data.frame(treatment, cellline, replicate) %>% 
#   mutate(cellline=as.factor(paste0(cellline, replicate))) %>% 
#   dplyr::select(-replicate)
# rownames(coldata)=colnames(df_numbers)
# 
# dds_fc <- DESeqDataSetFromMatrix(countData = df_numbers,
#                                  colData = coldata,
#                                  design = ~ cellline + treatment)
# 
# keep <- rowSums(counts(dds_fc) >= 10) >= 4
# dds_fc <- dds_fc[keep,]
# dds_fc$treatment <- relevel(dds_fc$treatment,ref="mock")
# group_combination <- c("treatment","CBE", "mock")

#####################
#####From Salmon#####
#####################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/christian/mouse_neurons_may_2020/salmon/")
files <- list.files(".", pattern = "*sf")
#files <- files[grepl('NSC', files)]
files <- files[grepl('BC', files)]
#files <- files[grepl('B6', files)]
#files <- files[!grepl('B62', files)]
strain <- rep(c("B6", "BC", "NSC"), each=6)
treatment <- files %>% strsplit("_") %>% sapply(.,'[', 2)
samples <- files %>% strsplit("_R1") %>% sapply(., '[', 1)
individual <-  files %>% strsplit("_") %>% sapply(.,'[', 1)
coldata <- data.frame(treatment, individual, row.names = samples)
coldata

###############
###Tx2gene ####
###############
tx_ids <- readr::read_tsv(file = files[[1]]) %>% .$Name 
tx_ids <- sapply(strsplit(tx_ids,"\\."), function(x) x[1])

tx2gene <- read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/mm10/mart_mm10_ensemblv100.txt") %>% 
  dplyr::select(`Transcript stable ID`, `Gene stable ID`)
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene, ignoreTxVersion = T)

# Multi Factor design to include sample information as a term
# This will account for differences between the samples while
# estimating the effect due to the treeatment
dds_salmon <- DESeqDataSetFromTximport(txi = txi_salmon,
                                          colData = coldata,
                                          design = ~ individual + treatment)

keep <- rowSums(counts(dds_salmon) >= 10) >= 4
dds_salmon <- dds_salmon[keep,]
dds_salmon$treatment <- relevel(dds_salmon$treatment,ref="mock")
group_combination <- c("treatment","CBE", "mock")

#########################
#####DATA EXPLORATION####
#########################
run_analysis(dds_salmon, group_combination, 1, 0.05, "mm10", explore_data = T)
#Dendogram
rld_salmon <- rlog(dds_salmon)
dists <- dist(t(assay(rld_salmon)))
plot(hclust(dists))

########################
### 2D PCA from DESeq ##
########################
#dds_norm_salmon <- DESeq(dds_salmon)
#norm_counts_salmon <- 1 + as.data.frame(counts(dds_norm_salmon, normalized=T))
DESeq2::plotPCA(rld_salmon, ntop=1000, intgroup = c("individual", "treatment"))
pca_data <- DESeq2::plotPCA(rld_salmon, intgroup = c("individual", "treatment"), returnData=T)
ggplot(pca_data,aes(x=PC1, y=PC2, color=name, shape=group)) +
  geom_point(size=5)

########################
##### Manual 2DPCA #####
########################
rld_salmon <- rld_salmon[, !(colnames(rld_salmon) %in% c("B62_mock", "B62_CBE", "BC1_mock", "BC1_CBE", "B61_mock", "B61_CBE"))]
rv <- rowVars(assay(rld_salmon))
select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pc <- prcomp(t(assay(rld_salmon)[select,]))
percentVar <- pc$sdev^2 / sum( pc$sdev^2 )

pca_matrix <- as.data.frame(pc$x)
pca_matrix <- as.data.frame(merge(pca_matrix, colData(dds_salmon), by="row.names", all.x=TRUE)) %>% column_to_rownames(var="Row.names")

ggplot(pca_matrix,aes(x=PC1, y=PC2, color=strain, shape=treatment)) +
  geom_point(size=7) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  theme(legend.title=element_blank(), text = element_text(size = 22), aspect.ratio=0.7)

######################
####### 3D PCA #######
######################
######################
library(rgl)
options(rgl.printRglwidget = F)
options(rgl.useNULL = F)
rgl.clear() 
rgl.open()

colors <- c("red", "blue")
plot3d(pc$x[,1:3], type='s', size = 3, col = colors, surface=FALSE, ellipsoid = TRUE)
text3d(pc$x[,1] + 2, pc$x[,2] + 2, pc$x[,3] + 2 , texts =rownames(pc$x), 
       color="black", family="serif",  font=5, cex=1)
#par3d(windowRect=c( 0,0,50,50 ))

setwd("~/Desktop/")

for (i in c(1:360)){
  rgl.viewpoint(i) 
}
rgl.viewpoint(300)
rgl.postscript("300.pdf", fmt = "pdf", drawText = T)


######################
###### dispersion ####
######################
#Plot of dispersion estimates over the average expression strength 
# First, gene-wise maximum likelihood estimates (MLE) are obtained using only the respective gene’s data (black dots). 
# Then, a curve (red) is fit to the MLEs to capture the overall trend of dispersion-mean dependence. 
# This fit is used as a prior mean for a second estimation round, which results in the final
# maximum a posteriori (MAP) estimates of dispersion (arrow heads). This can be understood as a shrinkage (along the blue arrows) 
# of the noisy gene-wise estimates toward the consensus represented by the red line. 
# The black points circled in blue are detected as dispersion outliers and not shrunk toward the prior (shrinkage would follow the dotted line).
# For clarity, only a subset of genes is shown, which is enriched for dispersion outliers. 
dds = DESeq(dds_salmon)
plotDispEsts(dds)

######################
###### pheatmap ######
######################
library("pheatmap")
#Select top30 most expressed genes over all samples
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:30]
df <- as.data.frame(colData(rld_salmon)[,c("strain","treatment")])
pheatmap(counts(dds, normalized=T)[select,], cluster_rows=F, show_rownames=F,
         cluster_cols=T, annotation_col=df, angle_col=45)

#Select top30 genes with largest variance 
rld_salmon <- rlog(dds_salmon)
select <- order(rowVars(assay(rld_salmon)), decreasing=T)[1:100]
pheatmap(assay(rld_salmon)[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=df, angle_col=45, fontsize = 15)

#Pairwise correlations
pheatmap(cor(assay(rld_salmon)), fontsize = 15, angle_col = 45)

###################
#### Outliers #####
###################
# Cooks's distance is a measure of how much a single sample is influencing
# the fitted coefficients for a gene, and a large value of cook's distance
# is intended to indicate an outlier count.
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), ylab="log10 (cook's distance)", range=0, las=2)

##################################################################
### Pairwise comparison of individual highest variance genes #####
##################################################################
pairwise_var <- function(individual){
  selected_ind <- assay(rld_salmon)[,grep(individual, colnames(assay(rld_salmon)))]
  select <- order(rowSds(selected_ind), decreasing = T)[1:100]
  out <- selected_ind[select,]
  return(rownames(out))
  
}
unique_individuals <- as.character(unique(rld_salmon$individual))
top_100k_pairwise <- as.data.frame(sapply(unique_individuals, pairwise_var))
top_100k_pairwise <- top_100k_pairwise[, -grep("B62", colnames(top_100k_pairwise))]
library(VennDiagram)
for (strain in c("B6", "BC", "NSC")){
  to_venn <- top_100k_pairwise[,grep(strain, colnames(top_100k_pairwise))]  
  grid.newpage()    
  if (ncol(to_venn) == 2){
    vectors <- list(to_venn[,1], to_venn[,2])
    alpha <- c(0.5, 0.5)
    fill <- c('lightyellow', 'paleturquoise')
    categories <- c("rep1", "rep3")
  } 
  else{
    vectors <- list(to_venn[,1], to_venn[,2], to_venn[,3])
    alpha <- c(0.5, 0.5, 0.5)
    fill <- c('lightyellow', 'paleturquoise', 'darkgrey')
    categories <- c("rep1", "rep2", "rep3")
  }
  plt <- venn.diagram(x = vectors,
                      na = "remove",
                      category.names = categories,
                      fill = fill,
                      alpha = alpha,
                      print.mode = c("raw"),
                      filename = NULL,
                      output = T,
                      main= strain
                      )
  grid::grid.draw(plt)
}

join <- function(x){
  joined <- left_join(as.data.frame(x), dplyr::select(ensembl_genes, c(gene_id, symbol, description)), 
                       by = c( "x" = "gene_id")) %>% distinct()
  show(head(joined))
  return(joined)
}
names(t)
t[[1]]
write_xlsx(t[[9]], path = paste0("NSC3.xlsx", sep=""), col_names=T, format_headers = T)
t <- apply(top_100k_pairwise, 2, join)
m <- do.call(rbind.data.frame, t)
out_sign_annot <- left_join(rownames_to_column(res), 
                            dplyr::select(ensembl_genes, c(gene_id, symbol, description)), 
                            by = c( "rowname" = "gene_id"))
######################
###### RUVSeq ########
######################
#Get empirical list of control genes to use in RUVseq in order to remove batch effects. 
#Select based on the lowest variance
library(RUVSeq)
empirical <- rownames(assay(rld_salmon)[rowVars(assay(rld_salmon)) < 0.1,])

colors = c(rep("red", times=5), rep("blue", times=21))
set <- newSeqExpressionSet(counts(dds_salmon))
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-2, 2))#, col=colors)
plotPCA(set, cex=1, col=colors)


par(mfrow = c(3, 2))
for(k in 1:6) {
  set_g <- RUVg(x = set, cIdx = empirical, k = k)
  plotPCA(set_g, cex=1.2, main = paste0('with RUVg, k = ',k))
}
normalized <- RUVg(set, empirical, k=3)
plot_umwanted_variation_factors(pData(normalized), groups=dds_salmon$individual, k=2)
plotRLE(set_g, outline=F, ylim=c(-2, 2))

#######################
######### DE ##########
#######################
log2cutoff <- 1
padjcutoff <- 0.05
setwd("~/Desktop/")

#################
#### salmon #####
#################
# This is a tricky example, where we have with grouped individuals (B6, BC and NSC),
# and we seek to test the group-specific effect of the CBE treatment while controlling 
# for individual effects (paired samples) 
dds_salmon$ind.n <- factor(rep(rep(1:3,each=2),3))

# Removing outlier sample (B62)
dds_salmon <- dds_salmon[, !(colnames(dds_salmon) %in% c("B62_mock", "B62_CBE"))]
dds_salmon$ind.n <- factor(c(as.character(rep(1:2,each=2)), as.character(rep(rep(1:3,each=2),2))))
design(dds_salmon) <- ~ strain + + strain:ind.n + strain:treatment
ml <- model.matrix( ~ strain + strain:ind.n + strain:treatment, colData(dds_salmon))
ml

# NSC
dds_NSC <- dds_salmon[,grep("NSC", colnames(dds_salmon))]
dds_NSC$individual <- droplevels(dds_NSC$individual)
# Removing NSC1 replicate for NSC analysis 
dds_NSC <- dds_salmon[,-grep("NSC1", colnames(dds_salmon))]

design(dds_NSC) <- ~ individual + treatment
ml <- model.matrix(~ individual + treatment, colData(dds_NSC))
ml

# B6
dds_B6 <- dds_salmon[, !(colnames(dds_salmon) %in% c("BC1_mock", "BC1_CBE", "NSC1_mock", "NSC1_CBE"))]
dds_B6$ind.n <- droplevels(dds_B6$ind.n)
ml <- model.matrix(~ strain + strain:ind.n + strain:treatment, colData(dds_B6))
ml
dds_B6 <- dds_salmon[, grep("B6", colnames(dds_salmon))]
dds_B6$individual <- droplevels(dds_B6$individual)
design(dds_B6) <- ~ individual + treatment
ml <- model.matrix(~ individual + treatment, colData(dds_B6))
ml

# BC
dds_BC <- dds_salmon[, !(colnames(dds_salmon) %in% c("BC1_mock", "BC1_CBE", "B61_mock",
                                                     "B61_CBE", "NSC1_mock", "NSC1_CBE"))]
dds_BC$ind.n <- droplevels(dds_BC$ind.n)
ml <- model.matrix(~ strain + strain:ind.n + strain:treatment, colData(dds_BC))
ml
dds_BC <- dds_salmon[, grep("BC", colnames(dds_salmon))]
dds_BC <- dds_BC[, -grep("BC1", colnames(dds_BC))]
dds_BC$individual <- droplevels(dds_BC$individual)
design(dds_BC) <- ~ individual + treatment
ml <- model.matrix(~ individual + treatment, colData(dds_BC))
ml[,-2]

# This design results in
ml <- model.matrix(~ strain + strain:ind.n + strain:treatment, colData(dds_salmon))
all.zero <- apply(ml, 2, function(x) all(x==0))
idx <- which(all.zero)
ml <- ml[,-idx]

ml

# Contrasts
group_combination <- list("strainNSC", "strainBC")

source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

list_de_salmon <- run_analysis(dds_salmon, group_combination,log2cutoff, padjcutoff, use_contrast = F, 
                               FC_shrinkage_method = "apeglm", 
                               genome = "mm10", annotate_locally = T, explore_data = F)

de_genes_salmon_apeglm_shrinkage <- rownames(list_de_salmon[[2]])
all_annot_salmon <- annotate_results(list_de_salmon[[1]], T, "mm10") 


as.data.frame(list_de_salmon[[1]])%>% rownames_to_column %>% filter(str_detect(rowname, 'ENSMUSG00000113149'))
as.data.frame(counts(dds_salmon), normalized=T) %>% rownames_to_column %>% filter(str_detect(rowname, 'ENSMUSG00000113149'))


###### Cross check if different approaches harbor the same genes ######
#With previous NOVOGENE analysis
novogene_DE_B6 <- readxl::read_excel("~/MEOCloud/analysis/christian/latest_sequecing_data/novogene_analysis/first_batch_mouse_neurons/results/4.DiffExprAnalysis/DEGlist/B6_CBE_vs_B6_Mock.genes.significant.DEA.xlsx") %>%
  pull(Gene)

B6_2samples_joined_analysis <- readxl::read_excel("~/Desktop/strainB6.treatmentCBE.xlsx") %>% pull(gene_id)
table(B6_2samples_joined_analysis %in% de_genes_salmon_apeglm_shrinkage)
BC_2samples_joined_analysis <-readxl::read_excel("~/Desktop/strainBC.treatmentCBE.xlsx") %>% pull(gene_id)
grid.newpage()    
plt <- venn.diagram(x = list(de_genes_salmon_apeglm_shrinkage, BC_2samples_joined_analysis),
                    #x = list(genes_old_analysis, de_genes_fc_apeglm_shrinkage, de_genes_salmon_apeglm_shrinkage)
                    na = "remove",
                    category.names = c("A", "B"),
                    fill = c('lightyellow', 'paleturquoise'),
                    alpha = c(0.5, 0.5),
                    #category.names = c("Old FC" , "FC normal shrinkage" , "Salmon normal shrinkage"),
                    #fill = c('lightyellow', 'paleturquoise', 'lightgreen'),
                    #alpha = c(0.5, 0.5, 0.5),
                    print.mode = c("raw"),
                    filename = NULL,
                    output = T)
grid::grid.draw(plt)


#####################################
######## Functional enrichment ######
#####################################
#GPROFILER
mock_cbe_gprofiler <- run_enrichment_gprofiler(de_genes_salmon_apeglm_shrinkage,custom_genes_background = rownames(all_annot_salmon),
                                               retrieve_only_sign_results = F, 
                                               exclude_iea = F, 
                                               retrieve_short_link = F,
                                               measure_under = F, 
                                               domain_scope = "annotated",
                                               organism = "mmusculus",
                                               sources= NULL)

head(mock_cbe_gprofiler %>% filter(source == "GO:MF"))
plot_enrichment_results(mock_cbe_gprofiler %>% dplyr::filter(significant == T), top_n = 10, short_term_size = F, add_info_to_labels = T, size_of_short_term = 500)

#FGSEA
test2 <- test %>% left_join(dplyr::select(m_df, human_gene_symbol, gene_symbol) %>% distinct(), by=c("symbol_mouse"="gene_symbol")) %>%
  test2 %>%
  dplyr::select(everything()) %>%
  summarise_all(funs(sum(!is.na(.))))


out <- run_fgsea_analysis(all_annot_salmon, 
                          #source="GO:BP",
                          is_ranked_already = F, 
                          gene_list_from = "DESeq2",
                          top_n = 30, 
                          npermutations = 5000,
                          organism = "mmusculus")

plot_enrichment_results(dplyr::filter(out, significant == T), 
                        rank_by = "NES", label = "Normalized Enrichment Score", top_n = 30, 
                        font_size = 12,
                        reverse_order = F,
                        add_info_to_labels =  T,
                        short_term_size = F, 
                        size_of_short_term = 50)


####################################################
###### DE cross check with previous analysis #######
####################################################
#With previous NOVOGENE analysis
novogene_DE_list <- readxl::read_excel("~/MEOCloud/analysis/christian/latest_sequecing_data/novogene_analysis/first_batch_mouse_neurons/results/4.DiffExprAnalysis/DEGlist/NSC_CBE_vs_NSC_Mock.genes.significant.DEA.xlsx") %>% 
  pull(Gene)
table(novogene_DE_list %in% de_genes_salmon_apeglm_shrinkage)

#With previous CBE analysis in fibroblasts
mart_mm10_hg38_map <- read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/mm10/mart_mm10_hg38_homologues.txt") %>% 
  dplyr::select(gene_id_human, gene_name_human, gene_id_mouse, gene_name_mouse, gene_description_human, gene_description_mouse) %>% distinct()

DE_mouse_human_map <- as_tibble(de_genes_salmon_apeglm_shrinkage) %>%
  rename(gene_id_mouse = value) %>%
  left_join(mart_mm10_hg38_map) %>% 
  drop_na()


fibroblasts_DE_list_fc <- readxl::read_excel("~/MEOCloud/analysis/christian/rna_seq_analysis/human_CBE_mock/diff_expression/DE_CBE_vs_mock_feature_counts.xlsx")                      
fibroblasts_DE_list_salmon <- readxl::read_excel("~/MEOCloud/analysis/christian/rna_seq_analysis/human_CBE_mock/diff_expression/DE_CBE_vs_mock_from_salmon.xlsx")                      

match <- DE_mouse_human_map[DE_mouse_human_map$gene_id_human %in% fibroblasts_DE_list_fc$gene_id,] 
DE_mouse_human_map[DE_mouse_human_map$gene_id_human %in% fibroblasts_DE_list_salmon$gene_id,]

DE_mouse_human_map[DE_mouse_human_map$gene_id_human %in% fibroblasts_DE_list_fc$gene_id,] %>% dplyr::select(gene_name_human, gene_description_human)

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
