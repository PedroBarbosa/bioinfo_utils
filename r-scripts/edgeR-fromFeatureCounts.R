source("http://bioconductor.org/biocLite.R")
library(edgeR)
setwd("C:/Users/User/Desktop")

##STEM##
stem_counts <- read.delim("C:/Users/User/Desktop/fromStringTie/edgeR-DESeq2/stem/featureCounts-out.txt", row.names = "Geneid",header=T)
stem_counts = stem_counts[,! colnames(stem_counts) %in% c("Geneid","Chr","Start","End","Strand","Length")]
##LEAF##
leaf_counts <- read.delim("C:/Users/User/Desktop/fromStringTie/edgeR-DESeq2/leaf/featureCounts-out.txt",row.names = "Geneid",header=T)
leaf_counts = leaf_counts[,! colnames(leaf_counts) %in% c("Geneid","Chr","Start","End","Strand","Length")]

counts <- list(stem_counts,leaf_counts)
groupletter <- c("c","f")
groupname <- c("stem", "leaf")
iter=1
for(letter in groupletter){

  s<-seq(1,ncol(stem_counts)/2)
  for (i in 1:length(s)){
    s[i]<-paste(letter,s[i],sep = "")
  }

  group=rep(s, each=2)
  factor(group)
  y <- DGEList(counts=counts[[iter]],group=group)
  
  
   #normalization
  keep <- rowSums(cpm(y)>1) >= 2
  cat("Number of genes keeped after low counts filtration:\n")
  cat(length(keep[keep==TRUE]),"\n")
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")
  y <- estimateDisp(y)
  
  
  #DE
  comb<-combn(unique(s),2)
  for (i in 1:ncol(comb)){
    combination <-paste(comb[,i], collapse = "-")
    print(paste("Processing ",combination, " pair.."))
    et <-exactTest(y, pair=comb[,i])
    topgenes <- topTags(et, n=Inf,sort.by = "p.value", p.value = 0.05)
    topgenes_FDR <- topgenes$table$FDR <= 0.05
    topgenes_filtered=topgenes[topgenes_FDR,]
    de <- summary(decideTestsDGE(et, p=0.05, adjust="BH"))
    colnames(de) <-paste("#",c(combination),sep="")
    rownames(de) <-c("#underexpressed","#no diffExp", "#overexpressed")
    outfile_down = paste("fromStringTie/edgeR-DESeq2/", combination, "_", groupname[iter],"_edgeR_down", ".txt", sep="")
    outfile_up = paste("fromStringTie/edgeR-DESeq2/", combination, "_", groupname[iter],"_edgeR_up", ".txt", sep="")
    write.table(de,outfile_down, sep = "\t", col.names=T,append=T, quote=F)
    write.table(de,outfile_up, sep = "\t", col.names=T,append=T, quote=F)
    if(nrow(topgenes_filtered)!=0) {
      col_vector=c("gene",colnames(topgenes_filtered))
      col_vector=c("gene",colnames(topgenes_filtered))
      columns=as.matrix(col_vector)
      t_columns=t(columns)
      t_columns[1,1] = paste("#",t_columns[1,1], sep="")
      
      ##underregulated
      under_expressed=subset(topgenes_filtered$table, logFC < 0)
      
      if(nrow(under_expressed) > 20000){
        cat("#Top under expressed genes:\n", file=outfile_down, append=T)
        write.table(t_columns,outfile_down,sep="\t",row.names= F,col.names = F, append=T,quote=F)
        write.table(head(under_expressed,20000),outfile_down, sep = "\t", col.names=F,append=T,quote=F)
      }
      else if (nrow(under_expressed) != 0){
        cat("#Top under expressed genes:\n", file=outfile_down, append=T)
        write.table(t_columns,outfile_down,sep="\t",row.names= F,col.names = F, append=T,quote=F)
        write.table(under_expressed,outfile_down, sep = "\t", col.names=F,append=T,quote=F)
      }
      else{
        cat("#No under expressed genes for this comparison!\n", file=outfile_down,append=T)
      }
      
     
      
      
      #upregulated
      over_expressed=subset(topgenes_filtered$table, logFC > 0)
      if(nrow(over_expressed) > 20000){
        cat("#Top over expressed genes:\n", file=outfile_up, append=T)
        write.table(t_columns,outfile_up,sep="\t",row.names= F,col.names = F, append=T,quote=F)
        write.table(head(over_expressed,20000),outfile_up, sep = "\t", col.names=F,append=T,quote,F)
      }
      else if (nrow(over_expressed) != 0){
        cat("#Top under expressed genes:\n", file=outfile_up, append=T)
        write.table(t_columns,outfile_up,sep="\t",row.names= F,col.names = F, append=T,quote=F)
        write.table(over_expressed,outfile_up, sep = "\t", col.names=F,append=T,quote=F)
      }
      else{
        cat("#No over expressed genes for this comparison!\n", file=outfile_up,append=T)
      }
      
    
    }
    else {
      cat("No genes with differential expression for this comparison!\n")
    }
    print("DONE!!")
  }
  iter<- iter + 1
}


































##VERIFICATION###
s<-seq(1,ncol(stem_counts)/2)
for (i in 1:length(s)){
  s[i]<-paste("c",s[i],sep = "")
}
group=rep(s, each=2)
factor(group)
y <- DGEList(counts=stem_counts,group=group)


#normalization
keep <- rowSums(cpm(y)>1) >= 2
cat("Number of genes keeped after low counts filtrartion:")
length(keep[keep==TRUE])

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y)


et <-exactTest(y, pair=c("c1","c3"))
topgenes <- topTags(et, n=Inf,sort.by = "p.value", p.value = 0.05)
topgenes_FDR <- topgenes$table$FDR <= 0.05
topgenes_filtered=topgenes[topgenes_FDR,]
de <- summary(decideTestsDGE(et, p=0.05, adjust="BH"))
de
colnames(de) <-c(combination)
rownames(de) <-c("underexpressed","no diffExp", "overexpressed")
outfile = paste("edgeR-output/edgeR_", groupname[1],"_",combination, ".csv")
write.table(de,outfile, sep = ",", col.names=T,append=T)

