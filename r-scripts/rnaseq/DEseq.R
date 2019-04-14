source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(dplyr)
library(DESeq2)
library(biomaRt)
biocLite("AnnotationDbi")
biocLite("biomaRt")

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)


countData=read.csv("C:/Users/User/Desktop/fromStringTie/edgeR-DESeq2/stem/featureCounts-out.txt", quote="", sep="\t", row.names = 1) %>% 
  dplyr::select(-c(Chr,Start,End,Strand,Length)) %>% 
  as.matrix()



# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData) > 1,]

#Metadata
groupletter <- c("c")
groupname <- c("stem")

cond <-seq(1,ncol(countData)/2)

for (i in 1:length(cond)){
  cond[i]<-paste(groupletter[iter],cond[i],sep = "")
}

group=rep(cond, each=2)
metadata = data.frame(sample_id=colnames(countData), condition=group)
comb<-combn(unique(cond),2)
#Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,colData=metadata,design=~condition)
dds = DESeq(dds)

for (i in 1:ncol(comb)){
  combination <-paste(comb[,i], collapse = "-")
  print(paste("Processing ",combination, " pair.."))
  v <- strsplit(combination,"-")[[1]]
  res = results(dds, contrast=c("condition", v[[1]], v[[2]]))
  res = res[order(res$pvalue),]
}



res
summary(res)

#Gene enrichment tests:
#RNA-enrich
#htsint
#ToPASeq topology based pathway analysis of RNA-seq data
#goseq
#SeqGSEA