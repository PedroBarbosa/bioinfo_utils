source("http://bioconductor.org/biocLite.R")
library('ballgown')
setwd("C:/Users/User/Desktop")

data_directory = system.file('extdata', package='ballgown')
data_directory_stem = paste(data_directory, "/stringtie-stem", sep = "")
data_directory_leaf = paste(data_directory, "/stringtie-leaf", sep = "")


directories_vector = c(data_directory_stem,data_directory_leaf)
groupletter <- c("c","f")
groupname <- c("stem", "leaf")
iter=1


for(directory in directories_vector){
  numb_files = length(list.dirs(path = data_directory_stem,recursive = FALSE))
  cond <-seq(1,numb_files/2)
  for (i in 1:length(cond)){
    cond[i]<-paste(groupletter[iter],cond[i],sep = "")
  }
  
  comb<-combn(unique(cond),2)
  for (i in 1:ncol(comb)){
    pattern_bg <-paste(toupper(comb[,i]), collapse = "")
    final_pat=paste(toupper(groupletter[iter]),"[",gsub("[A-Z]","",pattern_bg), "]",sep="")
    outfile = paste(getwd(), "/", pattern_bg, "_", groupname[iter],"_ballgown", ".txt", sep="")
  
    print(paste("Processing ",pattern_bg, " pair.."))
    bg = ballgown(dataDir=directory, samplePattern=final_pat, meas='all')
    pData(bg)=data.frame(id=sampleNames(bg), group=rep(comb[,i], each=2))
    
    bgresults = stattest(gown=bg, meas="FPKM", covariate="group", feature="transcript", getFC=TRUE)
    diffExp=sum(bgresults$pval < 0.05, na.rm = TRUE) #& bgresults$qval < 0.05, na.rm = TRUE)
    print(paste("Number of tested transcripts:",nrow(indexes(bg)$t2g)))
    print(paste("Number of diffExp transcripts in ",pattern_bg, " combination:", diffExp))
    
    subset <- subset(bgresults, bgresults$pval<0.05, na.rm=TRUE)
    if(nrow(subset)!=0) {
      print("Changing transcript names..")
      subset[, which(colnames(subset)=="id" )] <- sapply(subset[,which(colnames(subset)=="id" )], as.character) #remove factor from transcript ID
      #For loop to change the transcript id on the data frame
      for (index in 1:nrow(subset)) { 
        row = subset[index, ]
        new_transcript = paste("TCONS_",paste(rep(0,8-nchar(as.character(row$id))),collapse=""),as.character(row$id),sep = "")
        subset[index,which(colnames(subset)=="id" )] = new_transcript
      }
      print("Sorting table by pvalue..")
      subset2 = subset[order(subset$pval),]
      #subset2 = subset[order(-subset$fc),] #ORDER BY FOLD CHANGES [Descending ordere]
      
      #Remove feature column and create vector of columns to write output file
      subset2$feature <- NULL
      columns = as.matrix(colnames(subset2))
      t_columns=t(columns)
      t_columns[1,1] = paste("#",t_columns[1,1], sep="")
      
      print("Writing to file..")
      if(nrow(subset2) > 20000){
        cat("#Expressed transcripts:\t", nrow(subset2),"\n", file=outfile, append=T)
        write.table(t_columns,outfile,sep="\t", row.names= F,col.names = F,append=T, quote=F)
        write.table(head(subset2,20000),outfile, sep = "\t", row.names=F,col.names=F,append=T, quote=F)
      }
      else if (nrow(over_expressed) != 0){
        cat("#Expressed transcripts:\t", nrow(subset2),"\n", file=outfile, append=T)
        write.table(t_columns,outfile,sep="\t",row.names= F,col.names = F, append=T, quote=F)
        write.table(subset2,outfile, sep = "\t", row.names= F,col.names=F,append=T, quote=F)
      }
    }
    else {
      cat("No genes with differential expression for this comparison!\n")
    }
    print("DONE!!")
#    gene_expression = gexpr(bg)
#   bgresults_genes <- stattest(gown=bg,  meas="FPKM",covariate="group", feature="gene", getFC = TRUE, na.rm =FALSE)
#    diffExp_genes=sum(bgresults_genes$pval < 0.05 & bgresults_genes$qval < 0.05, na.rm = TRUE)
#    print(paste("Number of diffExp genes in ",pattern_bg, " combination:", diffExp_genes))
    
  }
  iter <- iter + 1
}




#bgresults = stattest(gown=bg, meas="FPKM", covariate="group", feature="transcript")
#hist(bgresults$qval, main='Ballgown p-values: c1-c3 comparison', col="grey", xlab='p-values, highly expressed transcripts')










##################2 GROUP COMPARISON#######################
# C1 - C3
bg_c1c3 = ballgown(dataDir=data_directory_stem, samplePattern='C[13]', meas='all')
pData(bg_c1c3)=data.frame(id=sampleNames(bg_c1c3), group=rep(c("c1","c3"), each=2))

bgresults_c1c3 = stattest(gown=bg_c1c3, meas="FPKM", covariate="group", feature="transcript")#getFC=TRUE)
hist(bgresults_c1c3$pval, main='Ballgown p-values: c1-c3 comparison', col="grey", xlab='p-values, highly expressed transcripts')
diffExp_c1c3=sum(bgresults_c1c3$pval < 0.05, na.rm = TRUE)# & bgresults_c1c3$qval < 0.05)
diffExp_c1c3
subset <- subset(bgresults_c1c3, bgresults_c1c3$pval<0.05, na.rm=TRUE)
hist(subset$qval, main='Ballgown q-values: c1-c3 comparison', col="grey", xlab='Range of q-values for the more significant transcripts')

under_expressed=subset(subset, fc > 0)

subset[, which(colnames(subset)=="id" )] <- sapply(subset[,which(colnames(subset)=="id" )], as.character) #remove factor from transcript ID
#For loop to change the transcript id on the data frame
for (index in 1:nrow(subset)) { 
  row = subset[index, ]
  new_transcript = paste("TCONS_",paste(rep(0,8-nchar(as.character(row$id))),collapse=""),as.character(row$id),sep = "")
  subset[index,which(colnames(subset)=="id" )] = new_transcript
}
subset2 = subset[order(-subset$fc),]

write.table(bgresults_c1c3,"teste.csv", sep=",")
new_transcript
bgresults$qval = p.adjust(bgresults$pval, 'fdr')
sum(bgresults$qval < 0.05, na.rm=TRUE)
bgresults

#gene based analysis. The number of fragments per million mapped reads (FPM) for each transcript in a gene is
#calculated by multiplying its FPKM value by its length in kb. The gene's FPM is
#the sum of the transcript level FPM for all the transcripts for that gene. Gene-level FPKM is obtained
#by dividing the gene FPM by the gene's length.
gene_expression = gexpr(bg_c1c3)
head(gene_expression)
bgresults_genes <- stattest(gown=bg_c1c3,  meas="FPKM", feature="gene")






##############EDGER from BALLGOWN tables############################
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
exon_counts=eexpr(bg)
colnames(exon_counts)
group=rep(c("c1","c2","c3","c4","c5"), each=2)
factor(group)
y <- DGEList(counts=exon_counts,group=group)
#normalization
y <- calcNormFactors(y, method="TMM")
design <- model.matrix(~group)
#estimate dispersion
y <- estimateDisp(y,design)
#calculate DE
et <- exactTest(y, pair=c("c1","c3"))
et_FDR = topTags(et, n=Inf,sort.by = "p.value", p.value = 0.01)
et_FDR(et_FDR$table$PValue < 0.01 & et_FDR$table$FDR < 0.01)
subset_de=subset(et_FDR$table, subset=et_FDR$table$PValue < 0.01 & et_FDR$table$FDR< 0.01)
nrow(subset_de)
