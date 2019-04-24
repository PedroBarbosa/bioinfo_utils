#library(tidyverse)
#library(corrplot)
library(heatmap.plus)

plotheatmap <- function(x) {
  heatmap.plus(x,
            main = "Genes correlation", # heat map title
            #density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(3,5),     # widens margins around plot
            col=hmcol,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="col",     # only draw a row dendrogram
            Colv=TRUE, 
            Rowv = FALSE # turn off column clustering
            #,scale="row"
  )   
}
args = commandArgs(trailingOnly=TRUE)
features_counts_data=read.table(args[1], header=TRUE,sep="\t")

df_numbers = features_counts_data[ , -which(names(features_counts_data) %in% c("gene_name"))]
rownames(df_numbers) = df_numbers[,'Geneid']
df_numbers[,'Geneid'] <- NULL
#df_numberic = mutate_all(df_numbers, function(x) as.numeric(as.character(x)))
cor(df_numbers)
#hc <- hclust(as.dist(df_numbers), method = "ward.D2")
#pdf("hclust.pdf")
#plot(hc)
#dev.off()
print(length(df_numbers))

Sample5k = df_numbers[sample(nrow(df_numbers), 5000), ]
t_n <- t(Sample5k)
print(head(t_n))
pdf("clustering.pdf")
clusters <- hclust(dist(t(df_numbers)))
plot(clusters)
#plotheatmap(data.matrix(df_numbers))
dev.off()


