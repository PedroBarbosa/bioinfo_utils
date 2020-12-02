library(DESeq2)
library(RUVSeq)
library(tidyr)
library(reshape2)
library(corrplot)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")


###################################
######FROM FEATURE COUNTS #########
###################################
setwd("/Users/pbarbosa/Dropbox/imm/projects/cajal_body_rnaseq/analysis/de_analysis/")
read_counts <- read.table("featureCounts.txt", row.names = "Geneid", sep="\t",header = TRUE)
read_counts <- read_counts[,! colnames(read_counts) %in% c("Geneid","gene_name","Chr","Start","End","Strand","Length")]
colnames(read_counts) <- sub('X.home.pedro.barbosa.cb_data.star_alignments.trim_galored.','', colnames(read_counts))
colnames(read_counts) <- sub('_Aligned.sortedByCoord.out.bam','', colnames(read_counts))
colnames(read_counts) <- sub('myC_time','t', colnames(read_counts))

#Filter genes where at least five reads are present in at least 2 samples
filter <- apply(read_counts, 1, function(x) length(x[x>5])>=2)
filtered <- read_counts[filter,]

coldata <- data.frame(timepoints, row.names=colnames(filtered))
dds <- DESeqDataSetFromMatrix(countData = filtered,
                              colData = coldata,
                              design = ~ timepoints)

#######################
##### FROM SALMON #####
#######################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/mirek_CajalBodies/new_data_spike-ins/salmon/gencode_v33/")
timepoints <- as.factor(rep(c(0, 12, 24, 48), each=3)) 
replicates <- as.factor(rep(c(1, 2, 3), times=4))
myc_activation <- as.factor(rep(c(0,1), times=c(3,9)))
coldata <- data.frame(timepoints, replicates, myc_activation)
rownames(coldata) <- paste(paste0("t",coldata$timepoints), coldata$replicates, sep="_")
files <- list.files(path=".", pattern = "*sf")
names(files) <- rownames(coldata)


# Import gene level counts by generating counts from abundance
tx2gene <- readr::read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33_with_spikein.txt")
tx_ids <- readr::read_tsv(file = "myC_time0_rep1_quant.sf") %>% .$Name 
# tx_ids <- sapply(strsplit(tx_ids,"\\."), function(x) x[1])
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID version`, `Gene stable ID`) ,
                                by = c("value" = "Transcript stable ID version")) %>% drop_na()

# gene level counts 
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txIn = T)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ replicates + timepoints)

###################################################
########  RUVseq spike ins control  ###############
###################################################
genes <- rownames(dds)[grep("^ENS", rownames(dds))]
spikes <- rownames(dds)[grep("^ERCC", rownames(dds))]
set <- newSeqExpressionSet(counts(dds))

#Data exploration
colors <- brewer.pal(4, "Set2")
colors <- as.factor(rep(colors, each=3))

# Normalize to sample size
plotRLE(set, outline=FALSE, col=colors, ylim=c(-4, 4))
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, col=colors, ylim=c(-4, 4))


#Remove variation using spike-in control genes
#Estimate factors of unwanted variation
#The normalized counts are indeed simply the residuals from ordinary least squares regression of logY on W
#(with the offset if needed). The output from RUV is actually exponentiated (and optionally rounded) 
#to be interpretable as "normalized counts". Please, note that these are intended just for 
#visualization/exploration, as the RUV model has been tested only on supervised problems, 
#and works better when W is inserted in the model with the original counts (rather than modeling the normalized counts).

# If using k=1, the estimated factors of unwanted variation refer mostly
# to large differences found at timepoint 0, probably related with variation
# intrinic to the cells. On the other hand, using k=2 we see that the major
# source of unwanted variation arises at timepoint 12, which is the first 
# timepoint upon treatment and most gene expression outliers caused by mYC
# activation may reside 
normalized <- RUVg(set, spikes, k=2)
plot_umwanted_variation_factors(pData(normalized), groups=dds$timepoints, ksize=2, groups_to_color=dds$replicates)

#PCAs raw vs normalized
EDASeq::plotPCA(set, k=2, labels=T, col=as.character(colors), cex=1.5)
# Problem with this is that the first factor is included in the normcounts
# Couldn't find a way to use just the 2nd factor
EDASeq::plotPCA(normalized, k=2, labels=T, col=as.character(colors), cex=1.5)

#rerun DEseq with the new design to reestimante the parametes and results
dds$W2 <- normalized$W_2
dds$W1 <- normalized$W_1
design(dds) <- ~ W2 + replicates + timepoints


###############################
###### Data exploration #######
### without RUVg analysis #####
###############################
explore_data_based_on_transformed_variance(dds, "timepoints", 
                                           transform_method = "rlog",
                                           additional_column_to_shape = "replicates",
                                           batch_effect_to_remove=dds$W2)  


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
dds <- run_analysis(dds, 
                   group_combination="", 
                   log2cutoff=log2cutoff,
                   padjcutoff=padjcutoff,
                   use_contrast = F,
                   FC_shrinkage_method = FC_shrinkage_method,
                   stat_test = stat_test,
                   reduced_formula="timepoints",
                   annotate_locally = T,
                   return_just_dds = T)

res_12_0 <- results(dds, name = "timepoints_12_vs_0", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_12_0 <- lfcShrink(dds, coef = "timepoints_12_vs_0", type=FC_shrinkage_method, res= res_12_0)
res_12_0 <- annotate_results(res_12_0, T, "hg38") 
res_sign_12_0 <- subset(res_12_0, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)

res_24_0 <- results(dds, name = "timepoints_24_vs_0", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_24_0 <- lfcShrink(dds, coef = "timepoints_24_vs_0", type=FC_shrinkage_method, res= res_24_0)
res_24_0 <- annotate_results(res_24_0, T, "hg38") 
res_sign_24_0 <- subset(res_24_0, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)

res_48_0 <- results(dds, name = "timepoints_48_vs_0", alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
res_48_0 <- lfcShrink(dds, coef = "timepoints_48_vs_0", type=FC_shrinkage_method, res= res_48_0)
res_48_0 <- annotate_results(res_48_0, T, "hg38")
res_sign_48_0 <- subset(res_48_0, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)


#####################################
######## Functional enrichment ######
#####################################
#GPROFILER
# 12
gprofiler_12_0 <- run_enrichment_gprofiler(rownames(res_sign_12_0),
                                               custom_genes_background = rownames(res_12_0),
                                               retrieve_only_sign_results = T, 
                                               exclude_iea = F, 
                                               retrieve_short_link = F,
                                               measure_under = F, 
                                               domain_scope = "annotated",
                                               sources= NULL)

plot_enrichment_results(gprofiler_12_0, top_n = 20, short_term_size = F, add_info_to_labels = T, size_of_short_term = 500, font_size = 12)

fgsea_12_0 <- run_fgsea_analysis(res_12_0, is_ranked_already = F, top_n = 30, npermutations = 5000)
plot_enrichment_results(dplyr::filter(fgsea_12_0, significant == T), 
                        rank_by = "log10padj", label = "Log10 padj", top_n = 30, 
                        font_size = 12,
                        reverse_order = F,
                        add_info_to_labels =  T,
                        short_term_size = F, 
                        size_of_short_term = 50)

# 24
gprofiler_24_0 <- run_enrichment_gprofiler(rownames(res_sign_24_0),
                                           custom_genes_background = rownames(res_24_0),
                                           retrieve_only_sign_results = T, 
                                           exclude_iea = F, 
                                           retrieve_short_link = F,
                                           measure_under = F, 
                                           domain_scope = "annotated",
                                           sources= NULL)
dim(gprofiler_24_0)
plot_enrichment_results(gprofiler_12_0, top_n = 30, short_term_size = F,
                        add_info_to_labels = T, size_of_short_term = 500, font_size = 12)

fgsea_24_0 <- run_fgsea_analysis(res_24_0, is_ranked_already = F, top_n = 30, npermutations = 5000)
fgsea_24_0 <- fgsea_24_0 %>% filter(!str_detect(term_name, "INSULIN_LIKE_GROWTH_FACTOR_IGF"))
plot_enrichment_results(dplyr::filter(fgsea_24_0, significant == T), 
                        rank_by = "log10padj", label = "Log10 padj", top_n = 30, 
                        font_size = 12,
                        reverse_order = F,
                        add_info_to_labels =  T,
                        short_term_size = F, 
                        size_of_short_term = 50)
# 48
gprofiler_48_0 <- run_enrichment_gprofiler(rownames(res_sign_48_0),
                                           custom_genes_background = rownames(res_48_0),
                                           retrieve_only_sign_results = T, 
                                           exclude_iea = F, 
                                           retrieve_short_link = F,
                                           measure_under = F, 
                                           domain_scope = "annotated",
                                           sources= NULL)
dim(gprofiler_48_0)
plot_enrichment_results(gprofiler_48_0, 
                        rank_by = "log10padj", top_n = 30, short_term_size = F, size_of_short_term = 500, font_size = 12)

fgsea_48_0 <- run_fgsea_analysis(res_48_0, is_ranked_already = F, top_n = 30, npermutations = 5000)
fgsea_48_0 <- fgsea_48_0 %>% filter(!str_detect(term_name, "INSULIN_LIKE_GROWTH_FACTOR_IGF"))
plot_enrichment_results(dplyr::filter(fgsea_48_0, significant == T), 
                        rank_by = "log10padj", label = "Log10 padj", top_n = 30, 
                        font_size = 12,
                        reverse_order = F,
                        add_info_to_labels =  T,
                        short_term_size = F, 
                        size_of_short_term = 50)

# all
all <- run_enrichment_gprofiler(unique(c(rownames(res_sign_12_0), rownames(res_sign_24_0),rownames(res_sign_48_0))),
                         custom_genes_background = unique(c(rownames(res_12_0), rownames(res_24_0), rownames(res_48_0))),
                         retrieve_only_sign_results = T, 
                         exclude_iea = F, 
                         retrieve_short_link = F,
                         measure_under = F, 
                         domain_scope = "annotated",
                         sources= NULL)

plot_enrichment_results(all, 
                        rank_by = "log10padj", top_n = 30, short_term_size = F, font_size=12, size_of_short_term = 500)


####################
##### LRT test #####
####################
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
stat_test <- "LRT"
dds_lrt = run_analysis(dds, 
                       group_combination="", 
                       log2cutoff=log2cutoff,
                       padjcutoff=padjcutoff,
                       use_contrast = F,
                       FC_shrinkage_method = FC_shrinkage_method,
                       stat_test = stat_test,
                       reduced_formula="W2 + replicates",
                       annotate_locally = T,
                       return_just_dds = F)


dds_lrt <- DESeq(dds, test="LRT", reduced = ~ W2 + replicates)
res_lrt <- results(dds_lrt)
sig_res_LRT <- res_lrt %>%
  data.frame() %>%
  filter(padj < 0.00001) 
sig_res_LRT <- annotate_results(sig_res_LRT, T, "hg38") 

gprofiler_lrt <- run_enrichment_gprofiler(sig_res_LRT,
                                           custom_genes_background = rownames(res_lrt),
                                           retrieve_only_sign_results = T, 
                                           exclude_iea = F, 
                                           retrieve_short_link = F,
                                           measure_under = F, 
                                           domain_scope = "annotated",
                                           sources= NULL)
plot_enrichment_results(gprofiler_lrt, 
                        rank_by = "log10padj", top_n = 30, short_term_size = T, font_size=12, size_of_short_term = 500)


#######################################
############ snRNA genes ##############
#######################################
#Normalized counts
dds_normalized = DESeq(dds)
cts <- counts(dds_normalized)#, normalized = T)
rownames(cts) <- str_replace(rownames(cts),
                             pattern = ".[0-9]+$",
                             replacement = "")
cts <- as_tibble(cts, rownames = "gene_ids")

#Raw counts
raw_cts <- as_tibble(read_counts, rownames = "gene_ids")
raw_cts$gene_ids <- str_replace(raw_cts$gene_ids,
                                     pattern = ".[0-9]+$",
                                     replacement = "")

#snRNA
gene_names <- c("RNAU1-1", "RNU2-1", "RNU4-1", "RNU5A-1", "RNU5B-1", "RNU6-1", "RNU11", "RNU12")
gene_ids <- c("ENSG00000206652", "ENSG00000274585", "ENSG00000200795","ENSG00000199568", 
               "ENSG00000200156", "ENSG00000206625", "ENSG00000274978", "ENSG00000276027")
df_snRNA <- tibble("gene_names" = gene_names, "gene_ids" = gene_ids)

norm_cts_snRNA <- left_join(df_snRNA, cts)
raw_cts_snRNA <- left_join(df_snRNA, raw_cts)

#####################################
######### CAJAL BODY GENES ##########
#####################################
all_interesting <- read.table("~/MEOCloud/analysis/cajal_bodies/all_interesting_genes.txt") %>% pull(V1)
cb_genes <- read.table('~/MEOCloud/analysis/cajal_bodies/CB_genes.txt') %>% pull(V1)
histone_genes <- read.table("~/MEOCloud/analysis/cajal_bodies/histone_genes.txt") %>% pull(V1)

de_genes_12_0 <- res_sign_12_0$symbol

paste("Number of de genes 12 vs 0:", length(de_genes_12_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_12_0)))
paste("Number of cajal body genes in the whole gene list:", length(which(cb_genes %in% res_12_0$symbol)))
# cb_genes not in tables:
# cb_genes[which(!cb_genes %in% res_12_0$symbol)]: "FRG1P"  "NHP2L1" "SMN" 
cb_genes_present_12_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_12_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_12_0)))
paste("Number of histone genes in the whole gene list:", length(which(histone_genes %in% res_12_0$symbol)))
histone_genes_present_12_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_12_0)]))
# histone genes not in tables:
# histone_genes[which(!histone_genes %in% res_12_0$symbol)] "CPSF73" "TROVE2"

de_genes_24_0 <- res_sign_24_0$symbol
paste("Number of de genes 24 vs 0:", length(de_genes_24_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_24_0)))
cb_genes_present_24_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_24_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_24_0)))
histone_genes_present_24_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_24_0)]))


de_genes_48_0 <- res_sign_48_0$symbol
paste("Number of de genes 48 vs 0:", length(de_genes_48_0))
paste("Number of genes in the cajal body list:", length(cb_genes))
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% de_genes_48_0)))
cb_genes_present_48_0 <- as.vector(na.omit(cb_genes[which(cb_genes %in% de_genes_48_0)]))

paste("Number of genes in the histone list:", length(histone_genes))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% de_genes_48_0)))
histone_genes_present_48_0 <- as.vector(na.omit(histone_genes[which(histone_genes %in% de_genes_48_0)]))

# Check genes from LRT
paste("Number of cajal body genes in the DE list:", length(which(cb_genes %in% sig_res_LRT$symbol)))
paste("Number of histone list in the DE list:", length(which(histone_genes %in% sig_res_LRT$symbol)))

df_cb_genes  <- as.data.frame(t(bind_rows(lapply(list(cb_genes_present_12_0, 
                                                     cb_genes_present_24_0, 
                                                     cb_genes_present_48_0), 
                                                 as.data.frame.list)))) %>% 
  mutate_if(is.character , replace_na, replace = "-")
colnames(df_cb_genes) <- c("12h", "24h", "48h")


df_histone_genes <-  as.data.frame(t(bind_rows(lapply(list(histone_genes_present_12_0, 
                                                           histone_genes_present_24_0, 
                                                           histone_genes_present_48_0), 
                                                      as.data.frame.list)))) %>% 
  mutate_if(is.character , replace_na, replace = "-")

colnames(df_histone_genes) <- c("12h", "24h", "48h")
write_xlsx(df_cb_genes, path = "CB_genes_DE.xlsx", col_names=T, format_headers = T)
write_xlsx(df_histone_genes, path = "histone_genes_DE.xlsx", col_names=T, format_headers = T)


###################################
#### Cross check with Splicing ####
###################################
setwd("/Users/pbarbosa/Desktop/lobo/MCFONSECA-NFS/mcfonseca/shared/mirek_CajalBodies/new_data_spike-ins/splicing/tools_overlap/")
genes_with_splicing <- read_tsv(file ="tools_overlap_gene_id_fetched_from_symbol.tsv") %>% dplyr::select(1,2) %>% 
  rename(gene_name = `#gene_name`) %>% 
  dplyr::select(gene_id, gene_name) %>%
  filter (! duplicated(gene_id)) %>%
  drop_na() %>%
  column_to_rownames("gene_id") 

ALL_DE <- unique(c(rownames(res_sign_12_0), rownames(res_sign_24_0),rownames(res_sign_48_0)))
splicing_obtained_and_de <- intersect(rownames(genes_with_splicing), ALL_DE)

gprofiler_isec_splicing_de <- run_enrichment_gprofiler(splicing_obtained_and_de,
                                          custom_genes_background = NULL,
                                          retrieve_only_sign_results = T, 
                                          exclude_iea = F, 
                                          retrieve_short_link = F,
                                          measure_under = F, 
                                          domain_scope = "annotated",
                                          sources= NULL)
plot_enrichment_results(gprofiler_isec_splicing_de, 
                        rank_by = "log10padj", top_n = 30, short_term_size = F, font_size=12, size_of_short_term = 500)

# Check for CB and histone genes 
paste("Number of cajal body genes in the Splicing list:", length(which(cb_genes %in% genes_with_splicing$gene_name)))
paste("Number of histone list in the Splicing list:", length(which(histone_genes %in% genes_with_splicing$gene_name)))
