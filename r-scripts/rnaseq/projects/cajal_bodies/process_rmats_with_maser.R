library(maser)
library(rtracklayer)

########### CB genes ##################
cb_genes <- read.table("/Users/pbarbosa/analysis/cajal_bodies/all_interesting_genes.txt")

########################################
############# maser objects #############
########################################
setwd("~/Desktop/t0_t12_rmats/")
maser_t12_t0 <- maser(".", c("t0h", "t12h"), ftype = "JC")
maser_cov_filt_t12_t0 <- filterByCoverage(maser_t12_t0, avg_reads = 10)
maser_top_t12_t0 <- topEvents(maser_cov_filt_t12_t0, fdr = 0.05, deltaPSI = 0.2)
print(maser_top_t12_t0)

setwd("~/Desktop/t0_t24_rmats/")
maser_t24_t0 <- maser(".", c("t0h", "t24h"), ftype = "JC")
maser_cov_filt_t24_t0 <- filterByCoverage(maser_t24_t0, avg_reads = 10)
maser_top_t24_t0 <- topEvents(maser_cov_filt_t24_t0, fdr = 0.05, deltaPSI = 0.2)
print(maser_top_t24_t0)

setwd("~/Desktop/t0_t48_rmats/")
maser_t48_t0 <- maser(".", c("t0h", "t48h"), ftype = "JC")
maser_cov_filt_t48_t0 <- filterByCoverage(maser_t48_t0, avg_reads = 10)
maser_top_t48_t0 <- topEvents(maser_cov_filt_t48_t0, fdr = 0.05, deltaPSI = 0.2)
print(maser_top_t48_t0)

setwd("~/analysis/cajal_bodies/rmats")
cb_genes_with_diff_splicing_events <- character()
mylist <- list('t12_t0'=maser_top_t12_t0,'t24_t0'=maser_top_t24_t0,'t48_t0'=maser_top_t48_t0)

for (i in seq_along(mylist)){
  cat("\n")
  print(paste0("Analyzing ", names(mylist)[i], " comparison."))
  for (event in c("SE","A5SS","A3SS","RI","MXE")){
    top = summary(mylist[[i]], type = event)
    write.table(top, quote = F, sep= "\t", row.names = F, file=paste0(names(mylist)[i], "_sign_events_", event, ".tsv"))
    if (nrow(subset(top, geneSymbol %in% cb_genes$V1)) > 0){
      target <- subset(top, geneSymbol %in% cb_genes$V1) 
      print(paste0("Number of CB genes with any ", event, " diff event in", names(mylist)[i], ": ", nrow(target)))
      cb_genes_with_diff_splicing_events <- append(cb_genes_with_diff_splicing_events, target$geneSymbol)
   } 
    else{
      print(paste0("No CB gene has any ", event, " diff event in", names(mylist)[i]))
    }
  }
}


unique(cb_genes_with_diff_splicing_events)

maser::display(RBM41, "SE")
plotGenePSI(RBM41, type = "SE", show_replicates = TRUE)
volcano(myc_top, fdr = 0.05, deltaPSI = 0.5, type = "SE")
dotplot(myc_top)

########### Map splicig events to protein features ########
gtf_path <- "~/analysis/cajal_bodies/gencode.v31.annotation.gtf.gz"
ens_gtf <- rtracklayer::import.gff(gtf_path)

#Only 4 genes in t12 vs t0 comparison
cb_genes_with_diff_splicing_events
SMNDC1 <- geneEvents(maser_top_events, "SMNDC1")
plotGenePSI(SMNDC1, type = "SE", show_replicates = TRUE)
SMNDC1_mapped <- mapTranscriptsToEvents(SMNDC1, ens_gtf,ncores = 4)
head(annotation(SMNDC1_mapped, "SE"))
SMNDC1_annot <- mapProteinFeaturesToEvents(SMNDC1_mapped, c("Domain_and_Sites", "Topology"), by="category",ncores = 4)
plotUniprotKBFeatures(SMNDC1_annot, "SE", event_id = 29023, gtf = ens_gtf, 
                      features = c("domain", "chain"), show_transcripts = TRUE,ncores = 4)

NOP56 <- geneEvents(maser_top_events, "NOP56")
plotGenePSI(NOP56, type = "A5SS", show_replicates = TRUE)
NOP56_mapped <- mapTranscriptsToEvents(NOP56, ens_gtf, ncores = 4)
head(annotation(NOP56_mapped, "A5SS"))
NOP56_annot <- mapProteinFeaturesToEvents(NOP56_mapped, c("Domain_and_Sites", "Topology"), by="category", ncores = 4)
plotUniprotKBFeatures(NOP56_annot, "A5SS", event_id = 1538, gtf = ens_gtf, 
                      features = c("domain", "chain"), show_transcripts = TRUE, ncores = 4)

LARP7 <- geneEvents(maser_top_events, "LARP7")
plotGenePSI(LARP7, type = "A5SS", show_replicates = TRUE)
LARP7_mapped <- mapTranscriptsToEvents(LARP7, ens_gtf,ncores = 4)
head(annotation(LARP7_mapped, "A5SS"))
LARP7_annot <- mapProteinFeaturesToEvents(LARP7_mapped, c("Domain_and_Sites", "Topology"), by="category",ncores = 4)
plotUniprotKBFeatures(LARP7_annot, "A5SS", event_id = 4449, gtf = ens_gtf, 
                      features = c("domain", "chain"), show_transcripts = TRUE,ncores = 4)

NELFCD <- geneEvents(maser_top_events, "NELFCD")
plotGenePSI(NELFCD, type = "MXE", show_replicates = TRUE)
NELFCD_mapped <- mapTranscriptsToEvents(NELFCD, ens_gtf,ncores = 4)
head(annotation(NELFCD_mapped, "MXE"))
NELFCD_annot <- mapProteinFeaturesToEvents(NELFCD_mapped, c("Domain_and_Sites", "Topology"), by="category",ncores = 4)
plotUniprotKBFeatures(NELFCD_annot, "MXE", event_id = 980, gtf = ens_gtf, 
                      features = c("domain", "chain"), show_transcripts = TRUE,ncores = 4)

for(i in 1:nrow(cb_genes)){
  g <- as.character(cb_genes[i,])
  g_events <- geneEvents(maser_top_events, g)
  srsf6_mapped <- mapTranscriptsToEvents(g_events, ens_gtf)
}
maser_top_events
S <- geneEvents(maser_top_events, "SNRPF")
SNRPF


