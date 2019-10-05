library(tximport)
library(GenomicFeatures)
library(DRIMSeq)
library(rnaseqDTU)
library(RUVSeq)
library(dplyr)
library(stageR)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")

####################
#####INPUT DATA#####
####################
setwd("/Users/pbarbosa/analysis/cajal_bodies/salmon/")
timepoints <- as.factor(rep(c(0, 12, 24, 48), each=3)) 
replicates <- as.factor(rep(c("rep1", "rep2", "rep3"), times=4))
myc_activation <- as.factor(rep(c(0,1), times=c(3,9)))
samples <- data.frame(timepoints, replicates, myc_activation)
samples$sample_id <- paste(paste0("t",samples$timepoints), samples$replicates, sep="_")
files <- list.files(path=".", pattern = "*sf")
names(files) <- samples$sample_id


################################
#Import transcript level counts#
#by generating counts from abundance
################################
txi <- tximport(files, type = "salmon", txIn = T, txOut = T, 
                countsFromAbundance = "scaledTPM")

#cts <- txi$counts
#cts <- cts[rowSums(cts) > 0, ]

##################################
###Transcript to gene mapping#####
##################################
#gtf <- "../gencode.v31.annotation.with.spikein.ids.changed.gtf.gz"
#txdb.filename <- "../gencode.v31.annotation.with.spikein.ids.changed.sqlite"
#txdb <- makeTxDbFromGFF(gtf)
#saveDb(txdb, txdb.filename)

#txdb <- loadDb(txdb.filename)
#txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
#tab <- table(txdf$GENEID)
#txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
#txdf$GENEID <- substr(txdf$GENEID, 1, 15)

#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#genedesc <-  getBM(attributes=c('ensembl_gene_id_version','external_gene_name','description'), 
#                   filters = 'ensembl_gene_id_version', values = txdf$GENEID, mart =ensembl)
#write.table(genedesc, quote= FALSE, row.names = FALSE, sep="\t",file="../ensembl_gene_map.csv")
gene_map <- read.table("../mart_export.txt", quote="", sep="\t", header=T)



###################################################
########  RUVseq spike ins control   ##############
###################################################

#FromTximport calculates an offset matrix for the user,
#and uses the estimated counts directly)
#The countsFromAbundance approach instead 
#scales the counts up and down per gene to account for 
#any change in gene length due to differential isoform usage. 
#Because this event is not that common or severe -- globally -- 
#it doesn't end up scaling the counts that much for most genes.

dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = samples, 
                                design = ~timepoints)

#filter transcripts where at least five reads are present in at least 2 samples
filter <- apply(counts(dds), 1, function(x) length(x[x>5])>=2)
filtered <- counts(dds)[filter,]

tx <- rownames(filtered)[grep("^ENS", rownames(filtered))]
#by using the filtered counts, we loose 10 spike-in transcript due to their low expression
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
format(cor(filtered[spikes,],method = "pearson"), digits = 4)
filtered[spikes,]
set <- newSeqExpressionSet(filtered)

#Data exploration
colors <- brewer.pal(8, "Set2")
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(set, col=colors, cex=1.2)

#Remove variation using spike-in control genes
#Estimate factors of unwanted variation
normalized <- RUVg(set, spikes, k=2)
plot_umwanted_variation_factors(pData(normalized), dds$timepoints, 2)

#Plot normalized counts
plotRLE(normCounts(normalized), outline=FALSE, ylim=c(-4, 4), col=colors)
plotPCA(normCounts(normalized), col=colors, cex=1.2)

#rerun DEseq with the new design to reestimante the parametes and results
samples$W1 <- normalized$W_1
samples$W2 <- normalized$W_2
#design(dds) <- ~ W2 + W1 + timepoints
#design(dds) <- ~ W1 + timepoints

############################
#######DRIMSeq analysis#####
############################
#Number of million of reads in each experiment
cts <- txi$counts
colSums(cts)/1e6
filter <- apply(cts, 1, function(x) length(x[x>5])>=2)
cts <- cts[filter,]

###########################################
####DRIMseq obj from count-scaled TPMs ####
###########################################
cts_filt <- cts[! grepl("ERCC", rownames(cts)),]
all(rownames(cts_filt) %in% gene_map$Transcript.stable.ID.version)
gene_map <- gene_map[match(rownames(cts_filt),gene_map$Transcript.stable.ID.version),]
counts <- data.frame(gene_id=gene_map$Gene.stable.ID.version,
                     feature_id=gene_map$Transcript.stable.ID.version,
                     cts_filt)
d <- dmDSdata(counts=counts, samples=samples)
d <- dmFilter(d,
              min_samps_feature_expr=3, min_feature_expr=10,
              min_samps_feature_prop=3, min_feature_prop=0.1,
              min_samps_gene_expr=12, min_gene_expr=10)

#####################################################
###DRIMseq obj with normalized counts from RUVseq####
#####################################################
cts_filt <- data.frame(normCounts(normalized))[! grepl("ERCC", rownames(normCounts(normalized))),]
gene_map <- gene_map[match(rownames(cts_filt),gene_map$Transcript.stable.ID.version),]

counts <- data.frame(gene_id=gene_map$Gene.stable.ID.version,
                     feature_id=gene_map$Transcript.stable.ID.version,
                     cts_filt)
d_with_normalized <- dmDSdata(counts=counts, samples=samples)
d_with_normalized <- dmFilter(d_with_normalized,
              min_samps_feature_expr=3, min_feature_expr=10,
              min_samps_feature_prop=3, min_feature_prop=0.1,
              min_samps_gene_expr=12, min_gene_expr=10)

#Number of genes with given number of filtered isoforms 
table(table(counts(d_with_normalized)$gene_id))


#######################################
########Create design matrix###########
#######################################
design_full_original <- model.matrix(~ timepoints, data=DRIMSeq::samples(d))
#design_full <- model.matrix(~ W1 + timepoints, data=DRIMSeq::samples(d))

#Apparently to model for any kind of change comparing to time0, but matrix has some formulation error (not full rank)
design_full <- model.matrix(~ myc_activation + timepoints, data = DRIMSeq::samples(d))

#To allow any pairwise comparison I want, then play with contrast argument.
design_full <- model.matrix(~ -1 + timepoints, data = DRIMSeq::samples(d))
design_full

d_with_normalized <- dmPrecision(d_with_normalized, design = design_full_original)
d_original_with_baseline <- dmPrecision(d, design= design_full_original)
d_reduced <- dmPrecision(d, design=design_full)

d_with_normalized <- dmFit(d_with_normalized, design=design_full_original, verbose=1)
d_original_with_baseline <- dmFit(d_original_with_baseline, design=design_full_original, verbose=1)
d_reduced <- dmFit(d_reduced, design=design_full, verbose=1)
design_full


register(MulticoreParam((7)))
head(coefficients(d), level = "feature")
#d <- dmTest(d, coef="timepoints0", verbose=1)

d <- dmTest(d_reduced, contrast=c(-1, 1, 0, 0), verbose=1)
d <- dmTest(d_with_normalized, coef="timepoints12", verbose=1)
d12_0 <- dmTest(d_with_normalized, coef="timepoints12", verbose=1)
d24_0 <- dmTest(d_with_normalized, coef="timepoints24", verbose=1)
d48_0 <- dmTest(d_with_normalized, coef="timepoints48", verbose=1)
samples(d)

################
####Results#####
################
get_results <- function(d, gene_map, threshold=0.05){
  no.na <- function(x) ifelse(is.na(x), 1, x)
  #Single p-value per gene, tests whether there is any 
  #differential transcript usage within the gene
  plotPValues(d, level="gene")
  res_per_gene <- DRIMSeq::results(d)
  res_per_gene <- res_per_gene[order(res_per_gene$pvalue, decreasing = FALSE), ]
  res_per_gene$pvalue <- no.na(res_per_gene$pvalue)
  res_per_gene <- inner_join(res_per_gene, unique(gene_map[, c("Gene.stable.ID.version","Gene.name",
                                                                    "Gene.description")]), by= c("gene_id" = "Gene.stable.ID.version"))
  res_per_gene_sign <- subset(res_per_gene, pvalue < threshold)
  
  #Single p-value per transcript, which tests whether the
  #proportions for this transcript changed within the gene
  plotPValues(d, level="feature")
  res_per_tx <- DRIMSeq::results(d, level="feature")
  res_per_tx <- res_per_tx[order(res_per_tx$pvalue, decreasing = FALSE), ]
  res_per_tx$pvalue <- no.na(res_per_tx$pvalue)
  res_per_tx <- inner_join(res_per_tx, unique(gene_map[,c("Gene.stable.ID.version","Gene.name",
                                                               "Gene.description")]), by= c("gene_id" = "Gene.stable.ID.version"))
  res_per_tx_sign <- subset(res_per_tx, pvalue < threshold)
  return(list(res_per_gene, res_per_tx, res_per_gene_sign, res_per_tx_sign))
}


############################
########stageR##############
############################
run_stageR <- function(drimseq_res, gene_map){
  res_per_gene <- drimseq_res[[1]]
  res_per_tx <- drimseq_res[[2]]
  
  pScreen <- res_per_gene$pvalue
  names(pScreen) <- res_per_gene$gene_id
  pConfirmation <- matrix(res_per_tx$pvalue, ncol=1)
  rownames(pConfirmation) <- res_per_tx$feature_id
  tx2gene <- res_per_tx[,c("feature_id", "gene_id")]
  #for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=FALSE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                    onlySignificantGenes=TRUE)
  })
  drim.padj <- subset(drim.padj, transcript < 0.05)
  drim.padj  <- drim.padj [order(drim.padj$transcript, decreasing = FALSE), ]  
  
#  inner_join(res_per_gene, unique(gene_map[, c("Gene.stable.ID.version","Gene.name",
#                                               "Gene.description")]), by= c("gene_id" = "Gene.stable.ID.version"))
  
  drim.padj <- inner_join(drim.padj, unique(res_per_tx[,c("gene_id", "Gene.name", "Gene.description")]), 
                          by = c("geneID" = "gene_id"))
  drim.padj <- drim.padj %>% rename(pvalue_adj_gene = gene, pvalue_adj_transcript = transcript)
  return (drim.padj)
}

########12_0#######
setwd("/Users/pbarbosa/analysis/cajal_bodies/salmon/")
t12_0 <- get_results(d12_0, gene_map)
drim.padj_12_0 <- run_stageR(t12_0, gene_map)
write.table(drim.padj_12_0, quote= FALSE, row.names = FALSE, sep="\t", file="DRIMseq_final_t12_0.tsv")

############24_0############
t24_0 <- get_results(d24_0, gene_map)
drim.padj_24_0 <- run_stageR(t24_0, gene_map)
write.table(drim.padj_24_0, quote= FALSE, row.names = FALSE, sep="\t", file="DRIMseq_final_t24_0.tsv")

###########48_0############
t48_0 <- get_results(d48_0, gene_map)
drim.padj_48_0 <- run_stageR(t48_0, gene_map)
write.table(drim.padj_48_0, quote= FALSE, row.names = FALSE, sep="\t", file="DRIMseq_final_t48_0.tsv")


####### MERGE #########
no.na.first <- function(x, y, z) ifelse(is.na(x), y, paste0(y, z))
no.na.second <- function(x,z){
  if (is.na(x[4])){
    return(x[2])
  }
  else if (is.na(x[2])){
    return(z)
  }
  else
    return(paste0(x[2],z))
}
m12_24 <- merge(drim.padj_12_0[,c("txID", "pvalue_adj_gene", "pvalue_adj_transcript")],
              unique(drim.padj_24_0[,c("txID", "pvalue_adj_gene", "pvalue_adj_transcript")]),
              by = "txID", all=T)

m12_24$timepoint = ""
m12_24$timepoint <- no.na.first(m12_24$pvalue_adj_gene.x, m12_24$timepoint, "t12")
m12_24$timepoint <- no.na.first(m12_24$pvalue_adj_gene.y, m12_24$timepoint, "t24")

m12_24_48 <- merge(m12_24[,c("txID", "timepoint")],
                unique(drim.padj_48_0[,c("txID", "pvalue_adj_gene", "pvalue_adj_transcript")]),
                by = "txID", all=T)


m12_24_48$timepoint <- apply(m12_24_48, 1, no.na.second, z= "t48")
m12_24_48 <- inner_join(m12_24_48[,c("txID", "timepoint")], unique(gene_map[,c("Transcript.stable.ID.version","Gene.stable.ID.version",
                                                      "Gene.name", "Gene.description")]),
                        by= c("txID" = "Transcript.stable.ID.version"))


m12_24_48 <- m12_24_48[order(m12_24_48$timepoint, m12_24_48$Gene.stable.ID.version, decreasing = FALSE), ]
m12_24_48 <- m12_24_48[,c(1,3,4,5,2)]
write.table(m12_24_48, quote= FALSE, row.names = FALSE, sep="\t", file="DRIMseq_final_merged.tsv")


#####################################
#######CROSS CHECK WITH DE GENES#####
#####################################
cb_genes <- read.table("/Users/pbarbosa/analysis/cajal_bodies/all_interesting_genes.txt")
cb_formation_genes_de <- read.table("/Users/pbarbosa/analysis/cajal_bodies/CB_genes_DE.csv", header = T, sep="\t") 
histone_genes_de <- read.table("/Users/pbarbosa/analysis/cajal_bodies/histone_genes_DE.csv", header = T, sep="\t")

cb_vect <- c(t(cb_formation_genes_de))
cb_uniq <- unique(cb_vect[cb_vect != ""])
hist_vect <- c(t(histone_genes_de))
hist_uniq <- unique(hist_vect[hist_vect != ""])
all_cb_de_genes <- unique(c(hist_uniq, cb_uniq))

de_genes_all_t12_t0 <- read.table("/Users/pbarbosa/analysis/cajal_bodies/T12_vs_T0_all_annotated.csv", sep="\t")
de_genes_all_t24_t0 <- read.table("/Users/pbarbosa/analysis/cajal_bodies/T24_vs_T0_all_annotated.csv", sep="\t")
de_genes_all_t48_t0 <- read.table("/Users/pbarbosa/analysis/cajal_bodies/T48_vs_T0_all_annotated.csv", sep="\t")

cb_genes_with_any_dtu <- character()
cb_de_genes_with_any_dtu <- character()
cb_noDE_genes_but_with_any_dtu <- character()
cb_DE_withNO_DTU <- character()
for(i in 1:nrow(cb_genes)){
  g <- as.character(cb_genes[i,])
  if(g %in% m12_24_48$Gene.name){
    cb_genes_with_any_dtu <- append(cb_genes_with_any_dtu, g)
    #print(res_per_gene_sign[grep(g,res_per_gene_sign$external_gene_name),])
    if(g %in% all_cb_de_genes){ #genes with DTU and DGE
      cb_de_genes_with_any_dtu <- append(cb_de_genes_with_any_dtu, g)
    }
    else{ #genes with DTU but no DGE
      cb_noDE_genes_but_with_any_dtu <- append(cb_noDE_genes_but_with_any_dtu, g)
    }
  }
  else if(g %in% all_cb_de_genes){ #genes with DE but no DTU
    print(g)
    cb_DE_withNO_DTU <- append(cb_DE_withNO_DTU,g)
  }
}


nrow(cb_genes)
length(cb_genes_with_any_dtu)
length(cb_de_genes_with_any_dtu)
length(cb_noDE_genes_but_with_any_dtu)
length(cb_DE_withNO_DTU)


plot_proportions <- function(gene, gene_map){
  pdf(paste0(gene,"_ribbon.pdf"))
    print(gene)
    gene_to_plot <- unique(gene_map[grep(paste0('^',gene,'$'), gene_map$Gene.name),]$Gene.stable.ID.version)
    print(as.character(gene_to_plot))
    plot <- plotProportions(d, gene_id = as.character(gene_to_plot), 
                  group_variable = "timepoints",
                  plot_type = "ribbonplot")
    plot(plot)
  dev.off()
  pdf(paste0(gene,"_bxpltWithDisp.pdf"))
  plot <- plotProportions(d, gene_id = as.character(gene_to_plot), 
                          group_variable = "timepoints",
                          plot_type = "boxplot2")
  plot(plot)
  dev.off()
}


# ENSG00000105197.11 example of t48 gene coef=timepoint48
# contrast=c(0, -1, 1, 1) , best pvalue ENSG00000150991.15 (no change)
plotProportions(d, gene_id = "ENSG00000153187.20", 
                group_variable = "timepoints",
                plot_type = "ribbonplot")

setwd("/Users/pbarbosa/analysis/cajal_bodies/salmon/plots/cb_genes_with_any_dtu")
for (g in as.character(cb_genes_with_any_dtu)){
  tryCatch( plot_proportions(g, gene_map), error = function(e) {print(paste("ERROR in gene", g))})
}

setwd("/Users/pbarbosa/analysis/cajal_bodies/salmon/plots/de_genes_with_no_dtu")
for (g in as.character(cb_DE_withNO_DTU)){
    tryCatch( plot_proportions(g, gene_map), error = function(e) {print(paste("ERROR in gene", g))})
}






  ##################################
###############DEXseq#############
##################################
library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)