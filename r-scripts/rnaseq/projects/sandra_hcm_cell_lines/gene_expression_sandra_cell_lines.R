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
length(isec)
write(isec, file = "~/Desktop/isec.txt")

