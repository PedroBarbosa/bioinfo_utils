library(maser)
args = commandArgs(trailingOnly=TRUE)
label1 = args[[1]]
label2 = args[[2]]
min_avg_reads = as.numeric(args[[3]])
fdr_threshold = as.numeric(args[[4]])
deltaPSI_threshold = as.numeric(args[[5]])

ftype = args[[6]] #JC or JCEC files
if (length(args) > 6) {
  gtf = args[[7]]
  plot_transcripts=TRUE 
  } else {
  plot_transcripts=FALSE 
}

if (ftype == "JC"){
  ftype_long = "JC.txt"
} else if (ftype == "JCEC"){
  ftype_long = "JCEC.txt"
}
#########################
####Coverage filter######
#########################
files <- list.files(path = ".", pattern=ftype_long)
get_mean <- function(x){
  inc_s1 <- base::sapply(strsplit(as.character(x["IJC_SAMPLE_1"]), ",", fixed=T), function(x) base::mean(as.numeric(x)))
  exc_s1 <- base::sapply(strsplit(as.character(x["SJC_SAMPLE_1"]), ",", fixed=T), function(x) base::mean(as.numeric(x)))
  inc_s2 <- base::sapply(strsplit(as.character(x["IJC_SAMPLE_2"]), ",", fixed=T), function(x) base::mean(as.numeric(x)))
  exc_s2 <- base::sapply(strsplit(as.character(x["SJC_SAMPLE_2"]), ",", fixed=T), function(x) base::mean(as.numeric(x)))
  return(min(c(inc_s1, exc_s1, inc_s2, exc_s2)))
}

for (f in files){
  inclusion_table <- read.table(f, header = T)
  nrow_before <- nrow(inclusion_table)
  show(paste0("Filtering ", f, " file.. (", nrow_before, " events)"))
  filtered_table <- inclusion_table[apply(inclusion_table, 1, get_mean) > min_avg_reads,]
  nrow_after <- nrow(filtered_table)
  show(paste0("Number of events kept by coverage filter (min reads ", min_avg_reads, "):", 
              nrow_after, "(", round(nrow_after/nrow_before*100,2), "%)"))
  if(!dir.exists("coverage_filt/")){
    dir.create("coverage_filt/")
  }
  write.table(filtered_table, file = paste0("coverage_filt/", f), sep = "\t", quote=F)
}


########################################
############# maser object #############
########################################
maser_obj <- maser("coverage_filt/", c(label1, label2), ftype = ftype)

####################
######Volcano#######
####################
as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
volcano_function <- function(type, maser_obj, fdr_threshold, deltaPSI_threshold){
  pdf(paste0("volcano_", type, ".pdf"))
  p <- maser::volcano(maser_obj, type = type, fdr = 0.05, deltaPSI = 0.2)
  plot(p)
  dev.off()
}
base::lapply(as_types, volcano_function, maser_obj, fdr_threshold, deltaPSI_threshold)

#############################
####dPSI/FDR filtering#######
#############################
maser_top_events <- maser::topEvents(maser_obj, fdr = fdr_threshold, deltaPSI = deltaPSI_threshold)
print(maser_top_events)

#############################
#####Plot transcripts########
#############################
if(plot_transcripts){
  ens_gtf <- rtracklayer::import.gff(gtf)
  dir.create("tx_plots")
  plot_transcript_tracks <- function(x, maser_obj, gtf_f, type) {
    ID=gsub(" ", "", x["ID"], fixed = TRUE)
    gene=x["geneSymbol"]
    gene_maser <- geneEvents(maser_obj, geneS = gene)
    pdf(paste0("tx_plots/",gene,"_",type,"_",ID,"_tx.pdf"))
    plotTranscripts(gene_maser, type = type, event_id = ID,
                    gtf = gtf_f, zoom = F, show_PSI = TRUE)
    dev.off()
  }
  
}


##########################
######Write output########
##########################
for (event in c("SE","A5SS","A3SS","RI","MXE")){
  tryCatch(
    {
      top = maser::summary(maser_top_events, type = event)
      write.table(top, quote = F, sep= "\t", row.names = F, file=paste0("sign_events_", event, ".tsv"))
      if(plot_transcripts){
        base::apply(top, 1, plot_transcript_tracks,maser_top_events,ens_gtf,event)  
      }
    },
    error=function(e){
      show(paste0("No ", event, " events passed the threshold"))
    }
  )
}