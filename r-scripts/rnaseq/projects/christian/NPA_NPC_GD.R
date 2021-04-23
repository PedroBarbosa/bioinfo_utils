library(DESeq2)
library(dplyr)
library(tidyr)
library(biomaRt)
library(RColorBrewer)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

#####################
#####From Salmon#####
#####################
setwd("/Users/pbarbosa/Desktop/NFS_Carmo/christian/human_fibroblasts/npa_npc_gd_april_2020/salmon")
files <- list.files(".", pattern = "*sf")

# Remove older CTRLs
# files <- files[!grepl("^S", files)]

# Remove bad sample
files <- files[!grepl("GM13205", files)]

files <- files[grepl("NPA|NPC|CTRL", files)]
files <- files[grepl("S|GD", files)]

names <- unlist(lapply(strsplit(files, split = ".sf"), `[`, 1))
groups <-  as.factor(unlist(lapply(strsplit(names, split = "_"), `[`, 2)))
coldata <- data.frame(groups)
coldata$batch <- as.factor(ifelse(grepl("GD|S", names), 1, 2))
coldata$phenotype <- as.factor(ifelse(grepl("CTRL", coldata$groups), "CTRL", "Disease"))
row.names(coldata) <- names

# Tximport
tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_ensemblv103.txt")
tx_ids <- readr::read_tsv(file = "GD1_GD.sf") %>% .$Name 
tx_ids <- sapply(strsplit(tx_ids,"\\."), function(x) x[1])
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID`, `Gene stable ID`) , by = c("value" = "Transcript stable ID")) %>% drop_na()
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene, ignoreTxVersion = T)

# Deseq2
dds <- DESeqDataSetFromTximport(txi = txi_salmon,
                                colData = coldata,
                                design = ~ groups)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds$groups <- relevel(dds$groups, ref="CTRL")


#########################
#####DATA EXPLORATION####
#########################
explore_data_based_on_transformed_variance(dds, "groups", "rlog", blind=T)
rld <- rlog(dds)

#Dendogram
dists <- dist(t(assay(rld)))
plot(hclust(dists))

########################
### 2D PCA from DESeq ##
########################
DESeq2::plotPCA(rld, ntop=500, intgroup = "groups")
pca_data <- DESeq2::plotPCA(rld, intgroup = "groups", returnData=T)

########################
##### Manual 2D PCA ####
########################
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))
percentVar <- pc$sdev^2 / sum( pc$sdev^2 )

pca_matrix <- as.data.frame(pc$x)
pca_matrix <- as.data.frame(merge(pca_matrix, colData(dds), by="row.names", all.x=TRUE))
pca_matrix$name <-sapply(strsplit(pca_matrix$Row.names, "_"), function(x) x[1])

ggplot(pca_matrix, aes(x=PC1, y=PC2, color=groups, shape=name)) +
  geom_point(size=5, stroke=1.5) + 
  scale_shape_manual(values=seq(1,15)) +
  xlab(paste0("PC1: ",round(percentVar[2] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[3] * 100),"% variance")) +
  coord_fixed() +
  theme(legend.title=element_blank(), text = element_text(size = 20))


###########################
####Correlation plot ######
###########################
##Multiple plot
# rld_salmon <- as.data.frame(assay(rld_salmon))
# combinations_mock <- as_tibble(combinations(ncol(rld_salmon)/2, 2, colnames(rld_salmon)[grepl("mock", colnames(rld_salmon))]))
# combination_CBE <- as_tibble(combinations(ncol(rld_salmon)/2, 2, colnames(rld_salmon)[grepl("CBE", colnames(rld_salmon))]))
# combinations <- rbind(combinations_mock, combination_CBE)
# 
# plotlist=list()
# i=1
# for(i in 1:nrow(combinations)) {
#   pair1 <- as.character(combinations[i, 1])
#   pair2 <- as.character(combinations[i, 2])
#   
#   #For RNA-seq counts, however, the expected variance grows with the mean. 
#   #For example, if one performs PCA directly on a matrix of counts
#   #or normalized counts (e.g. correcting for differences in sequencing depth), 
#   #the resulting plot typically depends mostly on the genes with highest counts because 
#   #they show the largest absolute differences between samples.
#   #A simple and often used strategy to avoid this is to take the logarithm of the normalized count values
#   #plus a pseudocount of 1; however, depending on the choice of pseudocount, 
#   #now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot,
#   #because taking the logarithm of small counts actually inflates their variance. 
#   #p<-ggplot(norm_counts_salmon, aes(x=log2(norm_counts_salmon[,pair1]), y=log2(norm_counts_salmon[, pair2]))) +
#   p<-ggplot(rld_salmon, aes(x=rld_salmon[,pair1], y=rld_salmon[, pair2])) +
#     geom_point(color='brown', alpha=0.7) +
#     labs(x=pair1, y=pair2) +
#     stat_cor(method="pearson", digits=3) + 
#     geom_smooth(method = "lm", colour="black", size=0.5) +
#     theme(text = element_text(size =15))
#   
#   plotlist[[i]] = ggplotGrob(p)
#   i= i + 1 
# } 
# grid.arrange(grobs=plotlist,nrow=2,ncol=3,
#              #bottom="Pairwise comparisons of gene expression (rlog transformation of count data)", 
#              gp=gpar(fontface="bold", col="YELLOW", fontsize=20))


##########################
####### DE analysis ######
##########################
log2cutoff <- 1.0
padjcutoff <- 0.05
FC_shrinkage_method <- "apeglm"
stat_test <- "Wald"
setwd('~/Desktop/')

#####################
##### Wald test #####
#####################
design(dds) <- ~ groups
design(dds) <- ~ batch + phenotype
#design(dds) <- ~ batch 
dds_ <- run_analysis(dds, 
                    group_combination="", 
                    log2cutoff=log2cutoff,
                    padjcutoff=padjcutoff,
                    use_contrast = F,
                    FC_shrinkage_method = FC_shrinkage_method,
                    stat_test = stat_test,
                    reduced_formula="phenotype",
                    annotate_locally = T,
                    return_just_dds = T)


res_gd_ctrl <- results(dds_, name = "groups_GD_vs_CTRL", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_gd_ctrl <- lfcShrink(dds_, coef = "groups_GD_vs_CTRL", type=FC_shrinkage_method, res= res_gd_ctrl)
res_gd_ctrl <- annotate_results(res_gd_ctrl, T, "hg38") 
res_sign_gd_ctrl <- subset(res_gd_ctrl, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
dim(res_sign_gd_ctrl)

res_npa_ctrl <- results(dds_, name = "groups_NPA_vs_CTRL", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_npa_ctrl <- lfcShrink(dds_, coef = "groups_NPA_vs_CTRL", type=FC_shrinkage_method, res= res_npa_ctrl)
res_npa_ctrl <- annotate_results(res_npa_ctrl, T, "hg38") 
res_sign_npa_ctrl <- subset(res_npa_ctrl, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
dim(res_sign_npa_ctrl)

res_npc_ctrl <- results(dds_, name = "groups_NPC_vs_CTRL", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_npc_ctrl <- lfcShrink(dds_, coef = "groups_NPC_vs_CTRL", type=FC_shrinkage_method, res= res_npc_ctrl)
res_npc_ctrl <- annotate_results(res_npc_ctrl, T, "hg38")
res_sign_npc_ctrl <- subset(res_npc_ctrl, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
dim(res_sign_npc_ctrl)


############################
## Results of single gene ##
############################
res_npa_ctrl %>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000117399'))
as.data.frame(counts(dds), normalized=T) %>% rownames_to_column %>% filter(str_detect(rowname, 'ENSG00000117399'))


##########################################
######## VENN PLOT between DE genes ######
##########################################
library(VennDiagram)

grid.newpage()    
plt <- venn.diagram(x = list(rownames(res_sign_gd_ctrl), rownames(res_sign_npa_ctrl), rownames(res_sign_npc_ctrl)),
                    na = "remove",
                    category.names = c("GD", "NPA", "NPC"),
                    fill = c('lightyellow', 'paleturquoise', 'lightgreen'),
                    alpha = c(0.5, 0.5, 0.5),
                    scaled=T,
                    cex = 1.5, cat.fontface = 2, lty =1,
                    print.mode = c("raw", "percent"),
                    filename = NULL,
                    output = T)
grid::grid.draw(plt)



#####################################
######## Functional enrichment ######
#####################################
#GPROFILER
res <- run_enrichment_gprofiler(rownames(res_sign_gd_ctrl),custom_genes_background = NULL,
                                               retrieve_only_sign_results = T, 
                                               exclude_iea = F, 
                                               retrieve_short_link = F,
                                               measure_under = F, 
                                               domain_scope = "annotated",
                                               sources= NULL)
dim(res)
plot_enrichment_results(res, top_n = 30, short_term_size = F,  font_size = 14, add_info_to_labels = T, size_of_short_term = 500)


#FGSEA
out <- run_fgsea_analysis(res_gd_ctrl, is_ranked_already = F, top_n = 40)
out_ <- out %>% dplyr::filter(source == "REAC")
plot_enrichment_results(dplyr::filter(out_, significant == T), 
                        rank_by = "NES", label = "Normalized Enriched Score", top_n = 40, 
                        font_size = 15,
                        reverse_order = T,
                        add_info_to_labels =  F,
                        short_term_size = F, 
                        size_of_short_term = 50)

plot_enrichment_results(dplyr::filter(out_, significant == T), 
                        rank_by = "NES", label = "Normalized Enriched Score",
                        top_n = 20, short_term_size = F, size_of_short_term = 500)



######################################################################
############## Heatmap of DE genes and overall PCA ###################
######################################################################
rld_de_genes <- assay(rld[rownames(res_sign_npc_ctrl),])
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
