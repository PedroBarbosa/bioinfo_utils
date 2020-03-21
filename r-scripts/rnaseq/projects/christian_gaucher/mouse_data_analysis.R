library(DESeq2)
library(tidyverse)
library(tidyr)
library(tximport)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

groups <- as.factor(c("CBE","CTRL", "CBE", "CTRL"))
replicates <- as.factor(c("B6","B6", "SV", "SV"))
samples <- data.frame(replicates, groups)
samples$sample_id <- paste(paste0(samples$replicates), samples$groups, sep="_")

###################################################
#################  FROM KALLISTO  #################
###################################################
setwd("/Users/pbarbosa/analysis/christian/standard_rna_seq_analysis/mouse/kallisto/")
files <- list.files(path=".", pattern = "*tsv")
names(files) <- samples$sample_id

tx2gene <- readr::read_tsv("B6_CBE_1.fq.gz_abundance.tsv") %>% 
  separate(target_id, into=c("tx_id", "gene_id_version"), sep = "\\|") %>% 
  dplyr::select(c(1,2)) %>%
  mutate(gene_id = str_replace(gene_id_version, pattern = ".[0-9]+$", replacement = "")) %>% 
  dplyr::select(1,3)
txi_kallisto <- tximport(files, type = "kallisto", txIn = TRUE, ignoreAfterBar = TRUE, tx2gene = tx2gene)

###################################################
#################  FROM SALMON  #################
###################################################
setwd("/Users/pbarbosa/analysis/christian/standard_rna_seq_analysis/mouse/salmon/")
files <- list.files(path=".", pattern = "*sf")
names(files) <- samples$sample_id

tx2gene <- readr::read_tsv("/Users/pbarbosa/analysis/genome_utilities/mm9/mart_mm9.txt")
tx_ids <- readr::read_tsv(file = "SV_CBE_1.fq.gz_quant.sf") %>% .$Name 
tx_ids <- tibble::enframe(name = NULL, x = tx_ids)
tx2gene <- tx_ids %>% left_join(dplyr::select(tx2gene, `Transcript stable ID version`, `Gene stable ID`) , by = c("value" = "Transcript stable ID version")) 
txi_salmon <- tximport(files, type = "salmon", txIn = T, tx2gene = tx2gene)



#########################################################
#################  FROM Feature counts  #################
########################################################
feature_counts <- read.table("/Users/pbarbosa/analysis/christian/standard_rna_seq_analysis/mouse/featureCounts_fromSTAR.txt",
                             header=T, sep="\t")

df_numbers = feature_counts[ , -which(names(feature_counts) %in% 
                                              c("Chr","Start", "End", "Strand", "Length", "gene_name"))]
colnames(df_numbers) <- gsub("X.home.pedro.barbosa.mcfonseca.shared.christian.mouse_stem_cells_sept_2018.star.",
                             "", colnames(df_numbers))
colnames(df_numbers) <- gsub("Aligned.sortedByCoord.out.bam",
                             "", colnames(df_numbers))
df_numbers = column_to_rownames(df_numbers, 'Geneid') 

coldata = data.frame(groups, row.names =  colnames(df_numbers))
deseq <- DESeqDataSetFromMatrix(countData = df_numbers,
                              colData = coldata,
                              design = ~ groups)

################################
#############DESeq2#############
################################
coldata = data.frame(groups, row.names =  colnames(txi_kallisto$counts))
deseq <- DESeqDataSetFromTximport(txi_salmon,
                                         colData = coldata,
                                         design = ~ groups)

keep <- rowSums(counts(deseq)) >= 10
deseq <- deseq[keep,]
dds = DESeq(deseq)
norm_counts = data.frame(counts(dds, normalized = T))

clusters <- hclust(dist(t(norm_counts)))
plot(clusters)
library(corrplot)
library(ggpubr)
library(gridExtra)
pairwise_comparisons <- vector(mode="list", length=6)
pairwise_comparisons[[1]] = c("B6_CBE", "B6_CTRL")
pairwise_comparisons[[2]] = c("B6_CBE", "SV_CTRL")
pairwise_comparisons[[3]] = c("B6_CBE", "SV_CBE")
pairwise_comparisons[[4]] = c("B6_CTRL", "SV_CTRL")
pairwise_comparisons[[5]] = c("B6_CTRL", "SV_CBE")
pairwise_comparisons[[6]] = c("SV_CTRL", "SV_CBE")
##Multiple plot
plotlist=list()
i=1
for (pair in pairwise_comparisons) {
  p<-ggplot(norm_counts, aes(x=norm_counts[,pair[1]], y=norm_counts[, pair[2]])) +
    geom_point(color='brown') +
    labs(x=pair[1], y=pair[2]) +
    stat_cor(method="pearson")
  
  plotlist[[i]] = ggplotGrob(p)
  i= i + 1 
}

grid.arrange(grobs=plotlist,nrow=2,ncol=3,bottom="Pairwise comparisons of the normalized gene expression values", 
             gp=gpar(fontface="bold", col="YELLOW", fontsize=15))


############ DE ###############
log2cutoff <- 1
padjcutoff <- 0.05
group_combination <- c("groups","CBE", "CTRL")
list_de <- run_analysis(dds, group_combination, log2cutoff, padjcutoff, "mm9", explore_data = T)


