library(DESeq2)
library(RUVSeq)
library(tidyverse)
library(tibble)
library(ggplot2)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

timepoints <- as.factor(rep(c(0, 12, 24, 48), each=3)) 
replicates <- as.factor(rep(c("rep1", "rep2", "rep3"), times=4))
samples <- data.frame(timepoints, replicates)
samples$sample_id <- paste(paste0("t",samples$timepoints), samples$replicates, sep="_")

###################################################
#################  FROM KALLISTO  #################
###################################################
setwd("/Users/pbarbosa/analysis/cajal_bodies/new_cb_data_with_spike_ins/kallisto/")
files <- list.files(path=".", pattern = "*tsv")
names(files) <- samples$sample_id

tx_ids <- readr::read_tsv("myC_time0_rep1_abundance.tsv") %>% 
  separate(target_id, into=c("tx_id", "gene_id_version"), sep = "\\|") %>% dplyr::select(c(1,2)) %>%
  mutate(gene_id = str_replace(gene_id_version, pattern = ".[0-9]+$", replacement = "")) %>% dplyr::select(1,3)

#Fix spike-in names 
tx2gene <- tx_ids %>% mutate(gene_id = coalesce(gene_id, tx_id))
txi_kallisto <- tximport(files, type = "kallisto", txIn = TRUE, ignoreAfterBar = TRUE, tx2gene = tx2gene)

###################################################
#################  FROM SALMON  ###################
###################################################
setwd("/Users/pbarbosa/analysis/cajal_bodies/new_cb_data_with_spike_ins/salmon/")
files <- list.files(path=".", pattern = "*sf")
names(files) <- samples$sample_id

tx2gene <- readr::read_tsv("/Users/pbarbosa/analysis/genome_utilities/hg38/mart_hg38.txt")
tx_ids <- readr::read_tsv(file = "myC_time0_rep1_quant.sf") %>% .$Name 
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID version`, `Gene stable ID`) , by = c("value" = "Transcript stable ID version")) 
#Fix spike-in names 
tx2gene <- tx2gene %>% mutate(gene_id = coalesce(`Gene stable ID`, `value`)) %>% dplyr::select(-`Gene stable ID`)
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene)


###################################################
###### Deseq2 dataset and RUVg normalization ######
###################################################
#KALLISTO
coldata = data.frame(timepoints, row.names =  colnames(txi_kallisto$counts))
dds_kallisto <- DESeqDataSetFromTximport(txi_kallisto,
                                colData = coldata,
                                design = ~ timepoints)

keep <- rowSums(counts(dds_kallisto)) >= 10
dds_kallisto <- dds_kallisto[keep,]

genes <- rownames(dds_kallisto)[grep("^ENS", rownames(dds_kallisto))]
spikes <- rownames(dds_kallisto)[grep("^ERCC", rownames(dds_kallisto))]
set_kallisto <- newSeqExpressionSet(counts(dds_kallisto))

colors <- rep(brewer.pal(4, "Set2"), each=3) 
plotRLE(set_kallisto, col=colors, outline=FALSE, ylim=c(-4, 4))
EDASeq::plotPCA(set_kallisto, col=colors, cex=1)
set_kallisto <- betweenLaneNormalization(set_kallisto, which="upper", offset = T)
plotRLE(set_kallisto, col=colors,outline=FALSE, ylim=c(-4, 4))
EDASeq::plotPCA(set_kallisto, col=colors, cex=1)

ruvg_normalized_kal <- RUVg(set_kallisto, spikes, k=3)
plot_umwanted_variation_factors(pData(ruvg_normalized_kal), dds_kallisto$timepoints, 3)
dds_kallisto$W3 <- ruvg_normalized_kal$W_3
design(dds_kallisto) <- ~ W3 + timepoints

#Remove batch effects using log2 transformed counts to explore data
#Did not work quite well.
dds_normalized_kallisto = DESeq(dds_kallisto)
vsd <- vst(dds_normalized_kallisto, blind=TRUE)
transformed_counts <- assay(vsd)
batch_removed <- limma::removeBatchEffect(transformed_counts, vsd$W4)
assay(vsd) <- batch_removed
DESeq2::plotPCA(vsd, intgroup=c("timepoints"))


###SALMON
coldata = data.frame(timepoints, row.names =  colnames(txi_salmon$counts))
dds_salmon <- DESeqDataSetFromTximport(txi_salmon,
                                         colData = coldata,
                                         design = ~ timepoints)

keep <- rowSums(counts(dds_salmon)) >= 10
dds_salmon <- dds_salmon[keep,]

genes <- rownames(dds_salmon)[grep("^ENS", rownames(dds_salmon))]
spikes <- rownames(dds_salmon)[grep("^ERCC", rownames(dds_salmon))]
set_salmon <- newSeqExpressionSet(counts(dds_salmon))

colors <- rep(brewer.pal(4, "Set2"), each=3) 
plotRLE(set_salmon, col=colors, outline=FALSE, ylim=c(-4, 4))
EDASeq::plotPCA(set_salmon, col=colors, cex=1)
set_salmon <- betweenLaneNormalization(set_salmon, which="upper", offset = T)
plotRLE(set_salmon, col=colors,outline=FALSE, ylim=c(-4, 4))
EDASeq::plotPCA(set_salmon, col=colors, cex=1)

ruvg_normalized_sal <- RUVg(set_salmon, spikes, k=4)
plot_umwanted_variation_factors(pData(ruvg_normalized_sal), dds_salmon$timepoints, 4)
dds_salmon$W3 <- ruvg_normalized_sal$W_3
design(dds_salmon) <- ~ W3 + timepoints
dds_normalized_salmon = DESeq(dds_salmon)

#######################################
##########snRNAs expression check######
#######################################
#It will be performed on the DESeq2 normalized gene counts,
#without removal of batch effects from myc induction. That, is taken into account when doing DE analysis
cts_kallisto <- counts(dds_normalized_kallisto, normalized = T)
cts_kallisto <- as_tibble(cts_kallisto, rownames = "gene_ids") 
cts_kallisto$t0 <- rowMeans(subset(cts_kallisto, select = c(2,3,4)), na.rm = FALSE)
cts_kallisto$t12 <- rowMeans(subset(cts_kallisto, select = c(5,6,7)), na.rm = FALSE)
cts_kallisto$t24 <- rowMeans(subset(cts_kallisto, select = c(8,9,10)), na.rm = FALSE)
cts_kallisto$t48 <- rowMeans(subset(cts_kallisto, select = c(11,12,13)), na.rm = FALSE)

cts_salmon <- counts(dds_normalized_salmon, normalized = T)
cts_salmon <- as_tibble(cts_salmon, rownames = "gene_ids") 
cts_salmon$t0 <- rowMeans(subset(cts_salmon, select = c(2,3,4)), na.rm = FALSE)
cts_salmon$t12 <- rowMeans(subset(cts_salmon, select = c(5,6,7)), na.rm = FALSE)
cts_salmon$t24 <- rowMeans(subset(cts_salmon, select = c(8,9,10)), na.rm = FALSE)
cts_salmon$t48 <- rowMeans(subset(cts_salmon, select = c(11,12,13)), na.rm = FALSE)


#ALL snRNA
snRNA_gene_names <- c("RNAU1-1", "RNU2-1", "RNU4-1", "RNU5A-1", "RNU5B-1", "RNU6-1", "RNU6-5", "RNU11", "RNU12")
snRNA_gene_ids <- c("ENSG00000206652", "ENSG00000274585", "ENSG00000200795","ENSG00000199568", 
              "ENSG00000200156", "ENSG00000206625", "ENSG00000206965", "ENSG00000274978", "ENSG00000276027")

#SNORD3A RNA
snRNA_gene_ids <- c("ENSG00000263934", "ENSG00000206652")
snRNA_gene_names <- c("SNORD3A", "RNAU1-1" )

#TOP 3 snRNA
snRNA_gene_names <- c("RNAU1-1", "RNU2-1", "RNU4-1")
snRNA_gene_ids <- c("ENSG00000206652", "ENSG00000274585", "ENSG00000200795")

#U2 pseudogene top mapped reads by STAR
snRNA_gene_names <- c("RNU2-1", "RNU2-2P", "RNU2-24P", "RNU2-22P", "RNU2-58P", "RNU2-65P", "RNU2-30P")
snRNA_gene_ids <- c("ENSG00000274585", "ENSG00000222328", "ENSG00000252639", "ENSG00000223198", "ENSG00000252212",
                    "ENSG00000222094", "ENSG00000252018")

#U2 pseudogene
snRNA_gene_names <- c("RNU2-1", "RNU2-2P")
snRNA_gene_ids <- c("ENSG00000274585", "ENSG00000222328")

df_snRNA <- tibble("gene_names" = snRNA_gene_names, "gene_ids" = snRNA_gene_ids)

#Averaged snRNAs expression over time
library(ggpubr)
library(RColorBrewer)

norm_cts_snRNA_kall <- left_join(df_snRNA, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(2:14)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_sal <- left_join(df_snRNA, cts_salmon) %>% drop_na() %>% dplyr::select(-c(2:14)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

p1_kallisto <- ggplot(norm_cts_snRNA_kall,  aes(x=timepoints, y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names), size=1) +
  geom_point(aes(color=gene_names), size=1.5) +
  scale_colour_brewer(palette="Set3")

p1_salmon <- ggplot(norm_cts_snRNA_sal,  aes(x=timepoints, y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names), size=1)+
  geom_point(aes(color=gene_names), size=1.5) +
  scale_colour_brewer(palette="Set3")

ggarrange(p1_kallisto, p1_salmon, ncol=2, nrow=1, common.legend = T, labels=c("Kallisto", "Salmon"), legend="bottom")


#Most expressed, check replicates expression
norm_cts_snRNA_with_replicates_kall <- left_join(df_snRNA, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(15:18))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_with_replicates_kall <- norm_cts_snRNA_with_replicates_kall %>% separate(timepoints, c("timepoint", "replicate"), "_")

norm_cts_snRNA_with_replicates_sal <- left_join(df_snRNA, cts_salmon) %>% drop_na() %>% dplyr::select(-c(15:18))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_with_replicates_sal <- norm_cts_snRNA_with_replicates_sal %>% separate(timepoints, c("timepoint", "replicate"), "_")

p2_kallisto <- ggplot(norm_cts_snRNA_with_replicates_kall, aes(x=timepoint, y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))

p2_salmon <- ggplot(norm_cts_snRNA_with_replicates_sal, aes(x=timepoint, y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))

ggarrange(p2_kallisto, p2_salmon, ncol=2, nrow=1, common.legend = TRUE, labels=c("Kallisto", "Salmon"), legend="bottom")

#################################
#######histone genes check ######
#################################
#HIST 1
hist_gene_names <- c("H1-1","H1-2","H1-3", "H1-4", "H1-5", "H1-6") 
hist_gene_ids <- c("ENSG00000124610", "ENSG00000187837","ENSG00000124575", "ENSG00000168298",
                   "ENSG0000018435", "ENSG00000187475")

#HIST 2
hist_gene_names <- c("H2AC18", "H2AC19","H2AC20", "H2AC21", "H3C13", "H3C14", "H3C15", "H3-2")
hist_gene_ids <- c("ENSG00000203812", "ENSG00000272196","ENSG00000184260","ENSG00000184270",
                  "ENSG00000183598", "ENSG00000203811", "ENSG00000203852", "ENSG00000273213")


df_histGenes <- tibble("gene_names" = hist_gene_names, "gene_ids" = hist_gene_ids)

norm_cts_histGenes_kall <- left_join(df_histGenes, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(3:14)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

norm_cts_histGenes_sal <- left_join(df_histGenes, cts_salmon) %>% drop_na() %>% dplyr::select(-c(3:14)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

p1_kallisto <- ggplot(norm_cts_histGenes_kall,  aes(x=timepoints, y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names), size=1)+
  geom_point(aes(color=gene_names), size=1.5) # +
  #scale_color_brewer(palette="Set3")

p1_salmon <- ggplot(norm_cts_histGenes_sal,  aes(x=timepoints, y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names), size=1)+
  geom_point(aes(color=gene_names), size=1.5) # +
  #scale_color_brewer(palette="Set3")


ggarrange(p1_kallisto, p1_salmon, ncol=2, nrow=1, common.legend = TRUE, labels=c("Kallisto", "Salmon"), legend="bottom")


#Two most expressed, check replicates expression
norm_cts_histGenes <- left_join(df_histGenes, cts) %>% drop_na() %>% dplyr::select(-c(15:18))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

norm_cts_snRNA_with_replicates <- norm_cts_histGenes %>% separate(timepoints, c("timepoint", "replicate"), "_")
p2 <- ggplot(norm_cts_snRNA_with_replicates, aes(x=timepoint, y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))
p2


#########################
#####DE ANALYSIS ########
#########################
log2cutoff <- 1
padjcutoff <- 0.05

#########################################
############# 12 vs 0 ###################
#########################################
group_combination <- c("timepoints","12", "0")
dds$timepoints <-relevel(dds$timepoints, ref="0")
list_de_12_0 <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "hg38", explore_data = FALSE)
rownames(list_de_12_0[[1]]) <- gsub("\\..*", "", rownames(list_de_12_0[[1]]))
all_annot_12_0 <- merge(as(list_de_12_0[[1]], "data.frame") , ensembl_genes, all.x = T, by=0)
#all_annot_12_0 <- annotate_results(list_de_8_0[[1]], "hg38") 

write.table(all_annot_12_0, quote= FALSE, row.names = TRUE, sep="\t",file="T12_vs_T0_all_annotated.csv")

fgseaResTidy <- run_fgsea_analysis(all_annot_12_0, "~/Downloads/c2.cp.reactome.v6.2.symbols.gmt", "Hallmark_pathways")
write.table(fgseaResTidy$pathway, quote = FALSE,row.names = FALSE,file="enriched_pathways_T8_vs_T0.csv")
fgseaResTidy <- run_fgsea_analysis(all_annot_12_0, "~/Downloads/c5.bp.v6.2.symbols.gmt", "GO_Biological_process")
goseq.results <- run_goseq_analysis(all_annot_12_0, "hg19", list_de_12_0[[2]]$Row.names)

