library(tximport)
mouse_mart <- read.csv("~/analysis/genome_utilities/mart_mm10_ensembl98.txt", header=T, sep="\t")
rownames(mouse_mart) <- mouse_mart$Transcript.stable.ID

setwd("~/analysis/christian/mouse")
files_mouse <- list.files(path=".", pattern = "*sf")
names(files_mouse) <- substr(files_mouse, start=1, stop=6)
txi_mouse <- tximport(files_mouse, type = "salmon", ignoreTxVersion=T, txIdCol=1, txOut = T, txIn = T)
rownames(txi_mouse$abundance) <- sub("\\..*","", rownames(txi_mouse$abundance))
head(rownames(txi_mouse$abundance))
head(rownames(mouse_mart))
merged_df <- merge(mouse_mart[,c("Gene.stable.ID", "Transcript.stable.ID",
                                 "Gene.name", "Gene.description")],
                   txi_mouse$abundance,
                   by=0)
write.table(merged_df, row.names = F, quote = F, sep = "\t", 
            file = "mouse_isoform_expression_target_genes.tsv")

###################
#####human########
human_mart <- read.csv("~/analysis/genome_utilities/mart_hg38.txt", header=T, sep="\t")
rownames(human_mart) <- human_mart$Transcript.stable.ID.version

setwd("~/analysis/christian/human")
files_human <- list.files(path=".", pattern = "*sf")
names(files_human) <- sub("\\.fq.*","",files_human)
txi_human<- tximport(files_human, type = "salmon", ignoreTxVersion=TRUE, txIdCol=1, txOut = T, txIn = T)
merged_df <- merge(human_mart[,c("Gene.stable.ID", "Transcript.stable.ID",
                                 "Gene.name", "Gene.description")],
                   txi_human$abundance,
                   by=0)

write.table(merged_df, row.names = F, quote = F, sep = "\t", 
            file = "human_isoform_expression_target_genes.tsv")
