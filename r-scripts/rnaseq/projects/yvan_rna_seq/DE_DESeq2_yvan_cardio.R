library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(RUVSeq)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")
ensembl_genes <- read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt") 
names(ensembl_genes)[names(ensembl_genes) == "Gene name"] <- "symbol"
names(ensembl_genes)[names(ensembl_genes) == "Gene description"] <- "description"
names(ensembl_genes)[names(ensembl_genes) == "Gene stable ID"] <- "gene_id"

cardio_genes <- read_tsv("/Users/pbarbosa/MEOCloud/analysis/yvan_rna_seq/cardiomyopathy_gene_panel_ids.txt", col_names = F)

############################
########INPUT DATA##########
############################
#FROM FEATURE COUNTS
setwd("/Users/pbarbosa/MEOCloud/analysis/yvan_rna_seq/diff_expression/deseq2_from_featurecounts/")
features_counts_data <- read.table("featureCounts.txt", 
                                header=TRUE,sep="\t")

cols_to_remove <- c("Chr", "Start", "End", "Strand", "Length", "gene_name")
str_to_remove_from_samples <- "X.home.pedro.barbosa.mcfonseca.shared.yvan_HCM_rnaseq.star_merged.uniq."
features_counts_data <- features_counts_data[, -which(names(features_counts_data) %in% cols_to_remove)]
                        
colnames(features_counts_data) <- gsub(pattern = str_to_remove_from_samples, "", colnames(features_counts_data)) 
colnames(features_counts_data) <- gsub(pattern = "_uniq.bam", "", colnames(features_counts_data))
features_counts_data <- column_to_rownames(features_counts_data, "Geneid")
#raw_counts <- features_counts_data[rowSums(features_counts_data)>50,]
keep <- rowSums(features_counts_data >= 10) >= 10
raw_counts <- features_counts_data[keep,]

########################
##### Manual 2DPCA #####
########################
groups <- as.factor(unlist(lapply(strsplit(colnames(raw_counts), split = "_"), `[`, 2)))
coldata <- data.frame(groups, row.names=colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = coldata,
                              design = ~ groups)
rld <- rlog(dds)
rv <- rowVars(assay(rld))

select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc <- prcomp(t(assay(rld)[select,]))
percentVar <- pc$sdev^2 / sum( pc$sdev^2 )

pca_data <- DESeq2::plotPCA(rld, intgroup = "groups", returnData=T)
pca_matrix <- as.data.frame(pc$x)
pca_matrix$groups <- pca_data$groups

colors <- c( "grey47 ","darkred", "darkblue") 
colors <- colors[as.numeric(pca_data$groups)]
ggplot(pca_matrix,aes(x=PC1, y=PC2, colour=groups)) +
  geom_point(size=5) + 
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  scale_colour_manual(values=setNames(colors, groups)) +
  theme(legend.title=element_blank(), text = element_text(size = 20))


###########################                                         
##########RUVSeq###########
###########################
groups <- as.factor(unlist(lapply(colnames(raw_counts), 
                                  function(x) {ifelse(grepl("Ctrl", x, fixed=F), "Ctrl", "Patient")})))

coldata <- data.frame(groups, row.names=colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = coldata,
                              design = ~ groups)

group_combination <- c("groups", "Patient", "Ctrl")

#Get empirical list of control genes to use in RUVseq in order to remove batch effects. Select based on the lowest p-value
list_de_sick_ctrl_temporary <- run_analysis(dds, group_combination = group_combination, 
                                            log2cutoff, padjcutoff, genome = "hg38", annotate_locally = T, explore_data = F)
top_5k_genes <- rownames_to_column(data.frame(list_de_sick_ctrl_temporary[[1]]), var="gene_id") %>%
  top_n(-5000, wt=padj) %>%
  arrange(pvalue) %>% pull(gene_id)
colors = c(rep("red", times=5), rep("blue", times=21))
set <- newSeqExpressionSet(counts(dds))
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-2, 2), col=colors)
plotPCA(set, cex=1, col=colors)

empirical <- rownames(set)[which(!(rownames(set) %in% top_2k_genes))]
hist(list_de_sick_ctrl_temporary[[1]][empirical,]$log2FoldChange)
par(mfrow = c(3, 2))
for(k in 1:6) {
  set_g <- RUVg(x = set, cIdx = empirical, k = k)
  plotPCA(set_g, col=colors, cex=1.2, main = paste0('with RUVg, k = ',k))
}
normalized <- RUVg(set, empirical, k=1)
plot_umwanted_variation_factors(pData(normalized),groups=dds$groups, k=4)
plotRLE(normalized, outline=F, ylim=c(-2, 2), col=colors)
plotRLE(set_g, outline=F, ylim=c(-2, 2), col=colors)
plotPCA(normalized, col=colors, cex=1)

#RUVseq doesn't seem to capture unwanted variation, as the PCAs after normalization mix CTRLs and Patients (different k)
#rerun DEseq with the new design to reestimante the parametes and results
#dds$W1 <- normalized$W_1
#design(dds) <- ~ W1 + groups


##########################
######Analysis############
##########################
log2cutoff <- 1
padjcutoff <- 0.05

#########################
#######Sick vs Ctrl #####
#########################
groups <- as.factor(unlist(lapply(colnames(raw_counts), 
                                               function(x) {ifelse(grepl("Ctrl", x, fixed=F), "Ctrl", "Patient")})))

coldata <- data.frame(groups, row.names=colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                       colData = coldata,
                                       design = ~ groups)

group_combination <- c("groups", "Patient", "Ctrl")

list_de_sick_ctrl <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, annotate_locally = T, genome = "hg38", explore_data = F)
#rownames(list_de_sick_ctrl[[1]]) <- gsub("\\..*$", "",rownames(list_de_sick_ctrl[[1]]))
#list_de_sick_ctrl[[1]] <- left_join(rownames_to_column(as.data.frame(list_de_sick_ctrl[[1]]), var = "gene_id"), 
#                                    dplyr::select(ensembl_genes, c("gene_id", "symbol"))) %>%
#  distinct() %>% column_to_rownames(var = "gene_id")

STRICT <- subset(list_de_sick_ctrl[[1]], padj<=padjcutoff & abs(log2FoldChange)>=2)


#########################
#######DCM vs Ctrl #####
#########################
groups <- as.factor(unlist(lapply(strsplit(colnames(raw_counts), split = "_"), `[`, 2)))
coldata <- data.frame(groups, row.names=colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                       colData = coldata,
                                       design = ~ groups)

group_combination <- c("groups", "DCM", "Ctrl")
list_de_dcm_ctrl <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, use_contrast = F, annotate_locally = T, genome = "hg38", explore_data = F)


#########################
#######ICM vs Ctrl ######
#########################
group_combination <- c("groups", "ICM", "Ctrl")
list_de_icm_ctrl <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, use_contrast = F, annotate_locally = T, genome = "hg38", explore_data = F)


#########################
#######DCM vs ICM #####
#########################
dds$groups <- relevel(dds$groups, ref = "DCM")
group_combination <- c("groups", "ICM", "DCM")
list_de_dcm_icm <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, use_contrast = F, annotate_locally = T, genome = "hg38", explore_data = F)
#No DE genes


#####################################
######## Functional enrichment ######
#####################################
#GPROFILER
de_genes <- rownames(list_de_sick_ctrl[[2]])
all_annot <- annotate_results(list_de_sick_ctrl[[1]], T, "hg38") 

sick_ctrl_gprofiler <- run_enrichment_gprofiler(de_genes,custom_genes_background = rownames(all_annot),
                                               retrieve_only_sign_results = F, 
                                               exclude_iea = F, 
                                               retrieve_short_link = F,
                                               measure_under = F, 
                                               domain_scope = "annotated",
                                               sources= NULL)
dim(sick_ctrl_gprofiler)
plot_enrichment_results(sick_ctrl_gprofiler, top_n = 20, short_term_size = F, add_info_to_labels = T, size_of_short_term = 500)

#GOseq
goseq.results <- run_goseq_analysis(de_genes, rownames(all_annot), "hg19")
plot_enrichment_results(dplyr::filter(goseq.results, significant == T), 
                        rank_by = "log10padj_over", top_n = 20, short_term_size = F, size_of_short_term = 200)

#FGSEA
out <- run_fgsea_analysis(all_annot, is_ranked_already = F, top_n = 20, npermutations = 2000)
plot_enrichment_results(dplyr::filter(out, significant == T), 
                        rank_by = "NES", label = "Normalized Enrichment Score", top_n = 20, 
                        font_size = 15,
                        reverse_order = T,
                        add_info_to_labels =  F,
                        short_term_size = F, 
                        size_of_short_term = 50)



########################################
###########CROSS-CHECK #################
########################################
de_genes_patient_ctrl <- list_de_sick_ctrl[[2]]$symbol
de_genes_dcm_ctrl <- list_de_dcm_ctrl[[2]]$symbol
de_genes_icm_ctrl <- list_de_icm_ctrl[[2]]$symbol

genes_splicing_heartDisease <-c("ACTG1","ANKS1A","NPPB","CACNA1C","CAMK2D","CD36","CD59","CLDND1","DNM2","ENG","ESRRG","FAM126A","FLNC","HDAC9",
"HOPX","IGF1","LDB3","LMO7","MRTFB","MYBPC3","MYH6","MYH7","PARVB","PDLIM3","RTN4","RYR2","SCN5A","TNNI3","TNNT2","TPM1","TRABD","TTN")

genes_splicing_heartDisease <- tibble::enframe(genes_splicing_heartDisease, name = NULL) %>% rename(gene_name = value) %>%
          left_join(genes_match, by=c("gene_name" = "Gene name")) %>%
          rename(gene_id = `Gene stable ID`) %>%
          dplyr::select(gene_id, gene_name) %>%
          filter (! duplicated(gene_id)) %>%
          column_to_rownames("gene_id") 

genes_heartDisease <- c("ABCC9","ACTC1","ACTG1","ACTN2","ALDH4A1","ANKRD1","ANKS1A","ATP5F1B","BAG3","BDH1","CAMK2D","CD36","CD59","COX17","CRYAB",
"CSRP3","CLDND1","CTF1","DCN","DES","CXCL17","DMD","DNAJC19","DSG2","DSP","EMD","EYA4","FHL2","FHOD3","FAM126A","FKTN","FLNC","GATAD1","GOT1"
,"HSP90AA1","HSP90AB1","HSPA5","ILK","LAMA4","LAMP2","LDB3","LMNA","LMO7","MYBPC3","MYH3","MYH6","MYH7","MYOZ1","MYOZ2","MYPN","NDUFV1",
"NEXN","PDHB","PDLIM3","PKP2","PLN","PSEN1","PSEN2","RBM20","RTN4","RYR2","SCN5A","SGCA","SGCD","SLC25A11","SLC25A13","SLC8A1","TAZ",
"TBX15","TCAP","TGFB1","TGFB3","TPM3","TMPO","TNNC1","TNNI3","TNNI3K","TNNT2","TPM1","TPM3","TPP1","TRDN-AS1","TTN","VCL"
)

genes_heartDisease <- tibble::enframe(genes_heartDisease, name = NULL) %>% rename(gene_name = value) %>%
  left_join(genes_match, by=c("gene_name" = "Gene name")) %>%
  rename(gene_id = `Gene stable ID`) %>%
  dplyr::select(gene_id, gene_name) %>%
  filter (! duplicated(gene_id)) %>%
  column_to_rownames("gene_id") 



cardio_genes_from_ncbi_panel <- cardio_genes %>%
  mutate(gene_id = gsub("\\..*","",.[[1]])) %>%
  dplyr::select(gene_id) %>% left_join(genes_match, by=c("gene_id" = "Gene stable ID")) %>%
  rename(gene_name = `Gene name`) %>%
  dplyr::select(gene_id,gene_name) %>%
  filter (! duplicated(gene_id)) %>%
  column_to_rownames("gene_id")

cardio_genes_all <- unique(c(genes_splicing_heartDisease$gene_name, genes_heartDisease$gene_name, cardio_genes_from_ncbi_panel$gene_name))


###################################
#####Gene with splicing events#####
###################################
setwd("/Users/pbarbosa/MEOCloud/analysis/yvan_rna_seq/")
genes_with_splicing <- read_tsv("overlap_counts_newest_vastools_appraoch.txt") %>% dplyr::select(1) %>% 
  rename(gene_name = `#gene`) %>% 
  left_join(genes_match, by=c("gene_name" = "Gene name")) %>%
  rename(gene_id = `Gene stable ID`) %>%
  dplyr::select(gene_id, gene_name) %>%
  filter (! duplicated(gene_id)) %>%
  drop_na() %>%
  column_to_rownames("gene_id") 

#duplicates
genes_with_splicing%>% rownames_to_column() %>% dplyr::filter(duplicated(gene_name))



splicing_obtained_and_de <- intersect(genes_with_splicing$gene_name, de_genes_patient_ctrl)
splicing_obtained_and_known_gene_splicing <- intersect(genes_with_splicing$gene_name,genes_splicing_heartDisease$gene_name )
cardio_genes_and_de <- intersect(cardio_genes_all, de_genes_patient_ctrl)

####################################
########DCM PAPER Heinig et.al #####
####################################
library(readxl)
#SPLICING GENES
dcm_large_paper <- read_excel("/Users/pbarbosa/Dropbox/imm/projects/yvan_cardioRNA-seq/patrícia/DCM_large_cohort_RNAseq_paper/DCM_vs_CTR_exon_usage_by_PSI_analysis.xlsx")
dcm_large_paper <- dcm_large_paper %>% dplyr::select(gene_id) %>% separate_rows(gene_id, sep="\\+") %>% pull(gene_id)
splicing_obtained_and_dcm_paper <- intersect(genes_with_splicing$gene_name, dcm_large_paper)
dplyr::select(dcm_large_paper, c(chrom,start,end,strand,gene_id,exonic_part_number)) %>% 
  arrange(chrom,start) %>% dplyr::select(-strand, strand) %>%
  write_tsv("/Users/pbarbosa/Dropbox/imm/projects/yvan_cardioRNA-seq/patrícia/DCM_large_cohort_RNAseq_paper/DiffSplicing_events.bed",
            col_names = F)
                             
#MYH6 MYH7 expression ratio
#MYH6 ratio against MYH7 is virtually absent in patients vs Ctrls (from paper)
ALL_GENES_DE <- list_de_sick_ctrl[[1]]
MYH6_MYH7 <- ALL_GENES_DE %>% filter(Gene_name == "MYH6" | Gene_name == "MYH7" )
#In our data, MYH6 is downregulated in patients, while MYH7 is not. Let's check now the expression on each individual and calculate the ratio
dds_norm <- DESeq(dds)
norm_expression <- counts(dds_norm, normalized=T)
target_gene_ids <- genes_match %>% filter(Gene_name == "MYH6" | Gene_name == "MYH7" ) %>% dplyr::select(Gene_stable_ID_version) %>% distinct() %>% pull()
to_plot <- as.data.frame((norm_expression[target_gene_ids,][1,] / (norm_expression[target_gene_ids,][1,] + norm_expression[target_gene_ids,][2,])))
colnames(to_plot)[1] <- "MYH6_MYH7_ratio"
ggplot(data = to_plot, mapping = aes(x = groups, y = MYH6_MYH7_ratio, fill=groups)) + 
  geom_boxplot() + 
  scale_y_continuous(name="ratio MYH6 / MYH7") 

#DE GENES
protein_coding <- read_excel("/Users/pbarbosa/Dropbox/imm/projects/yvan_cardioRNA-seq/patrícia/DCM_large_cohort_RNAseq_paper/DCM_vs_CTR_gene_expresison_ProteinCoding.xlsx") %>%
  mutate(gene_id =gsub("\\..*","",gene_id))
noncoding <- read_excel("/Users/pbarbosa/Dropbox/imm/projects/yvan_cardioRNA-seq/patrícia/DCM_large_cohort_RNAseq_paper/DCM_vs_CTR_gene_expresison_noncoding.xlsx") %>%
  mutate(gene_id =gsub("\\..*","",gene_id))

#Very few DE genes overlap.
DE_GENES_NO_NAs <- DE_GENES[!is.na(DE_GENES$symbol),]
library(VennDiagram)
venn.diagram( 
  x = list(rownames(STRICT), protein_coding$gene_id, noncoding$gene_id),
  category.names = c("DE" , "PC " , "NC"),
  filename = "DE_logFC2_venn_diagram_gene_IDs.png",
  output=TRUE
)


###RBM20 PAPER#####
rmb20_genes <- c("CAMK2D", "LDB3", "LMO7", "PDLIM3", "RTN4", "RYR2", "TTN")
splicing_obtained_and_rbm20_paper <- intersect(genes_with_splicing$gene_name, rmb20_genes)
splicing_obtained_and_rbm20_paper

RBM20_expression <- DESeq(dds)
op <- par(mar = c(8,4,2,2) + 0.1)
barplot(counts(RBM20_expression['ENSG00000203867.8'], normalized=T), las=2, col="grey")
par(op)


########################################
#######INTRON RETENTION ANALYSIS########
########################################
ALL_GENES_DE <- list_de_sick_ctrl[[1]]
ALL_GENES_DE <- ALL_GENES_DE[!is.na(ALL_GENES_DE$Gene_name),]
ALL_GENES_DE$rank_stat <- abs(ALL_GENES_DE$log2FoldChange / ALL_GENES_DE$lfcSE)
DE_GENES <- list_de_sick_ctrl[[2]]$symbol

MAJIQ <- read_tsv("../../majiq/deltapsi_analysis/deltaPSI0.2_p0.9_majiq_lsvs.csv") %>% drop_na(`#lsv_id`)
RMATS <- read_tsv("../../rmats/deltaPSI0.2_rmats.csv") %>% drop_na(`#coordinates_event_id`)
VASTOOLS <- read_tsv("../../vastools/yvan_dPSI_0.2_vastools.csv") %>% drop_na(`#Event_ID`)
PSICHOMICS <- read_tsv("../../psichomics/star_counts/input_data/all_concat_vs_Ctrl.csv")

majiq_RI_only_genes <- MAJIQ %>% filter(lsv_type == "RI") %>% dplyr::select(gene_name) %>% distinct() %>% pull
majiq_RI_only_otherEvents_genes <- MAJIQ %>% filter(str_detect(lsv_type, "RI") & !lsv_type == "RI") %>% dplyr::select(gene_name) %>% distinct() %>% pull
majiq_RI_with_otherEvents_genes <- MAJIQ %>% filter(str_detect(lsv_type, "RI")) %>% dplyr::select(gene_name) %>% distinct() %>% pull

rmats_RI_only_gene <- RMATS %>% filter(event_type == "RI") %>% dplyr::select(gene_name) %>% distinct() %>% pull()

vastools_RI_only_gene <- VASTOOLS %>% filter(str_detect(Major_type, "RI")) %>% dplyr::select(Gene) %>% distinct() %>% pull()

RI_genes <- as_tibble(unique(c(majiq_RI_with_otherEvents_genes, vastools_RI_only_gene,rmats_RI_only_gene)))

#many gene names from vastools disappeared
RI_genes <- RI_genes %>% left_join(., ALL_GENES_DE, by=c("value"="Gene_name")) %>% drop_na() %>% dplyr::select(c(value, rank_stat, log2FoldChange)) %>% rename("value"="Gene_name")
RI_genes$group <- "RI"
NO_RI_GENES <- anti_join(ALL_GENES_DE, RI_genes, by="Gene_name") %>% dplyr::select(c(Gene_name, rank_stat, log2FoldChange))
NO_RI_GENES$group <- "No_RI"

plot_df <- rbind(RI_genes, NO_RI_GENES)
plot_df$is_de <- plot_df$Gene_name %in% DE_GENES
plot_df$log2FoldChange = abs(plot_df$log2FoldChange)

ggplot(data = plot_df, mapping = aes(x = group, y = log2FoldChange, color = ifelse(is_de == T,"Yes",'No'))) + 
  geom_jitter(size= 0.9, alpha=0.7) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_x_discrete(labels=c("No intron retention","With intron retention"),name="Genes") +
  scale_y_continuous(name="Log2 Fold Change (absolute value)") +
  guides(color=guide_legend(title="Is differentially expressed")) +
  theme(text = element_text(size = 20)) +
  geom_boxplot(data = plot_df, mapping = aes(x = group, y = log2FoldChange, fill = group), alpha=0.6,inherit.aes = F) +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  guides(fill=F)



