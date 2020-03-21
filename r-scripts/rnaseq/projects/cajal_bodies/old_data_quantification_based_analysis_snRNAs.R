library(DESeq2)
library(RUVSeq)
library(tidyverse)
library(tibble)
library(ggplot2)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

timepoints <- as.factor(rep(c(0, 8, 24, 48, "wt"), each=2)) 
replicates <- as.factor(rep(c("rep1", "rep2"), times=5))
samples <- data.frame(timepoints, replicates)
samples$sample_id <- paste(paste0("t",samples$timepoints), samples$replicates, sep="_")

###################################################
#################  FROM KALLISTO  #################
###################################################
setwd("/Users/pbarbosa/analysis/cajal_bodies/old_cb_data/kallisto_old_data/")
files <- list.files(path=".", pattern = "*tsv")
names(files) <- samples$sample_id

tx_ids <- readr::read_tsv("0h_GH1742_R1.fq.gz_abundance.tsv") %>% 
  separate(target_id, into=c("tx_id", "gene_id_version"),sep = "\\|") %>% dplyr::select(c(1,2)) %>%
  mutate(gene_id = str_replace(gene_id_version, pattern = ".[0-9]+$", replacement = "")) %>% dplyr::select(1,3) %>%
  drop_na()

txi_kallisto <- tximport(files, type = "kallisto", txIn = TRUE, ignoreAfterBar = TRUE, tx2gene = tx_ids)

###################################################
#################  FROM SALMON  ###################
###################################################
setwd("/Users/pbarbosa/analysis/cajal_bodies/old_cb_data/salmon_old_data/")
files <- list.files(path=".", pattern = "*sf")
names(files) <- samples$sample_id

tx2gene <- readr::read_tsv("/Users/pbarbosa/analysis/genome_utilities/hg38/mart_hg38.txt")
tx_ids <- readr::read_tsv(file = "0h_GH1742_quant.sf") %>% .$Name 
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID version`, `Gene stable ID`) ,
                                by = c("value" = "Transcript stable ID version")) %>% drop_na()
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene)


###################################################
################# Deseq2 dataset ##################
###################################################
coldata = data.frame(timepoints, row.names =  colnames(txi_kallisto$counts))
dds_kallisto <- DESeqDataSetFromTximport(txi_kallisto,
                                colData = coldata,
                                design = ~ timepoints)

keep <- rowSums(counts(dds_kallisto)) >= 10
dds_kallisto <- dds_kallisto[keep,]
dds_kallisto_normalized = DESeq(dds_kallisto)
vsd_kallisto <- vst(dds_kallisto_normalized, blind=TRUE)
DESeq2::plotPCA(vsd_kallisto, intgroup=c("timepoints"))


dds_salmon <- DESeqDataSetFromTximport(txi_salmon,
                                         colData = coldata,
                                         design = ~ timepoints)
keep <- rowSums(counts(dds_salmon)) >= 10
dds_salmon <- dds_salmon[keep,]
dds_salmon_normalized = DESeq(dds_salmon)
vsd_salmon <- vst(dds_salmon_normalized, blind=TRUE)
DESeq2::plotPCA(vsd_salmon, intgroup=c("timepoints"))


#######################################
##########snRNAs expression check######
#######################################
#It will be performed on the DESeq2 normalized gene counts,
#without removal of batch effects from myc induction. That, is taken into account when doing DE analysis
cts_kallisto <- counts(dds_kallisto_normalized, normalized = T)
cts_kallisto <- as_tibble(cts_kallisto, rownames = "gene_ids") 
cts_kallisto$twt <- rowMeans(subset(cts_kallisto, select = c(2,3)), na.rm = FALSE)
cts_kallisto$t0 <- rowMeans(subset(cts_kallisto, select = c(4,5)), na.rm = FALSE)
cts_kallisto$t8 <- rowMeans(subset(cts_kallisto, select = c(6,7)), na.rm = FALSE)
cts_kallisto$t24 <- rowMeans(subset(cts_kallisto, select = c(8,9)), na.rm = FALSE)
cts_kallisto$t48 <- rowMeans(subset(cts_kallisto, select = c(10,11)), na.rm = FALSE)

cts_salmon<- counts(dds_salmon_normalized, normalized = T)
cts_salmon <- as_tibble(cts_salmon, rownames = "gene_ids") 
cts_salmon$twt <- rowMeans(subset(cts_salmon, select = c(2,3)), na.rm = FALSE)
cts_salmon$t0 <- rowMeans(subset(cts_salmon, select = c(4,5)), na.rm = FALSE)
cts_salmon$t8 <- rowMeans(subset(cts_salmon, select = c(6,7)), na.rm = FALSE)
cts_salmon$t24 <- rowMeans(subset(cts_salmon, select = c(8,9)), na.rm = FALSE)
cts_salmon$t48 <- rowMeans(subset(cts_salmon, select = c(10,11)), na.rm = FALSE)


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

norm_cts_snRNA_kall <- left_join(df_snRNA, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(2:12)) %>%
  pivot_longer(!starts_with("gene_names"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_sal <- left_join(df_snRNA, cts_salmon) %>% drop_na() %>% dplyr::select(-c(2:12)) %>%
  pivot_longer(!starts_with("gene_names"), values_to = "avg_norm_expression", names_to = "timepoints")

level_order <- c("twt", "t0", "t8", "t24", "t48")
p1_kallisto <- ggplot(norm_cts_snRNA_kall,  aes(x=factor(timepoints, level = level_order), y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names),size=1)+
  geom_point(aes(color=gene_names), size=1.5) +
  scale_colour_brewer(palette="Set3")

p1_salmon <- ggplot(norm_cts_snRNA_sal,  aes(x=factor(timepoints, level = level_order), y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names), size=1)+
  geom_point(aes(color=gene_names), size=1.5) +
  scale_colour_brewer(palette="Set3")

ggarrange(p1_kallisto, p1_salmon, ncol=2, nrow=1, common.legend = T, labels=c("Kallisto", "Salmon"), legend="bottom")



#Most expressed, check replicates expression
norm_cts_snRNA_with_replicates_kall <- left_join(df_snRNA, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(13:17))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_with_replicates_kall <- norm_cts_snRNA_with_replicates_kall %>% separate(timepoints, c("timepoint", "replicate"), "_")

norm_cts_snRNA_with_replicates_sal <- left_join(df_snRNA, cts_salmon) %>% drop_na() %>% dplyr::select(-c(13:17))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")
norm_cts_snRNA_with_replicates_sal <- norm_cts_snRNA_with_replicates_sal %>% separate(timepoints, c("timepoint", "replicate"), "_")

p2_kallisto <- ggplot(norm_cts_snRNA_with_replicates_kall, aes(x=factor(timepoint, level = level_order), y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))

p2_salmon <- ggplot(norm_cts_snRNA_with_replicates_sal, aes(x=factor(timepoint, level = level_order), y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))

ggarrange(p2_kallisto, p2_salmon, ncol=2, nrow=1, common.legend = T, labels=c("Kallisto", "Salmon"), legend="bottom")


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

norm_cts_histGenes_kall <- left_join(df_histGenes, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(2:12)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

norm_cts_histGenes_sal <- left_join(df_histGenes, cts_salmon) %>% drop_na() %>% dplyr::select(-c(2:12)) %>%
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

p1_kallisto <- ggplot(norm_cts_histGenes_kall,  aes(x=factor(timepoints, level=level_order), y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names),size=1)+
  geom_point(aes(color=gene_names),size=1.5)

p1_salmon <- ggplot(norm_cts_histGenes_sal,  aes(x=factor(timepoints, level=level_order), y=avg_norm_expression, group=gene_names)) +
  geom_line(aes(color=gene_names),size=1)+
  geom_point(aes(color=gene_names),size=1.5)

ggarrange(p1_kallisto, p1_salmon, ncol=2, nrow=1, common.legend = T, labels=c("Kallisto", "Salmon"), legend="bottom")

#Two most expressed, check replicates expression
norm_cts_histGenes <- left_join(df_histGenes, cts_kallisto) %>% drop_na() %>% dplyr::select(-c(13:17))  %>% 
  pivot_longer(starts_with("t"), values_to = "avg_norm_expression", names_to = "timepoints")

norm_cts_snRNA_with_replicates <- norm_cts_histGenes %>% separate(timepoints, c("timepoint", "replicate"), "_")
p2 <- ggplot(norm_cts_snRNA_with_replicates, aes(x=factor(timepoint, level=level_order), y=avg_norm_expression, fill=gene_names)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.75))
p2

