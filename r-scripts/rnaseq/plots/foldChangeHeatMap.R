library("readxl")
library("gplots")
library("RColorBrewer")
library("dplyr")

plotheatmap <- function(x) {
  heatmap.2(x,
            main = "Genes correlation", # heat map title
            #density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(3,5),     # widens margins around plot
            col=hmcol,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="col",     # only draw a row dendrogram
            Colv=TRUE, 
            Rowv = TRUE # turn off column clustering
            #,scale="row"
  )   
}

##FERRIER DATA##
data <- read_excel("TLX3_TALL_dif_gene_expression.xlsx")
gnames <- data[,"Gene"]
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- unlist(gnames, use.names=FALSE)
log2 <- log2(mat_data)
log2[log2==-Inf] <- -7

hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)
my_palette <- colorRampPalette(c("red", "yellow", "blue"))(n = 299)
plotheatmap(log2)


##TAMARA DATA##
data_NPA <- read_excel("DGE_C_NPA_vs_I_NPA.xlsx")
data_WT <- read_excel("DGE_C_WT_vs_C_NPA.xlsx")

gnames_NPA <- data_NPA[,c("Gene Name")]
npa_data <- data.matrix(data_NPA[,"DGE_C_NPA_vs_I_NPA_log2(FC)"])
colnames(npa_data) <- c("C_NPA_vs_I_NPA")
rownames(npa_data) <- unlist(gnames_NPA, use.names=FALSE)

gnames_WT <- data_WT[,c("Gene Name")]
wt_data <- data.matrix(data_WT[,"DGE_C_NPA_vs_I_NPA_log2(FC)"])
colnames(wt_data) <- c("C_WT_vs_C_NPA")
rownames(wt_data) <- unlist(gnames_WT,use.names=FALSE)
                     
merged <- merge(wt_data, npa_data, by = 0, all = TRUE)
gnames <- merged[,1]
final <- data.matrix(merged[,2:ncol(merged)])
rownames(final) <- gnames
final[is.na(final)] <- 0


df <- as.data.frame(final)
lowest_WT <- df[order(df$C_WT_vs_C_NPA),][0:15,]
highest_WT <- df[order(df$C_WT_vs_C_NPA, decreasing = TRUE),][0:15,]

lowest_NPA <- df[order(df$C_NPA_vs_I_NPA),][0:15,]
highest_NPA <- df[order(df$C_NPA_vs_I_NPA, decreasing = TRUE),][0:15,]

top_FC <- rbind(lowest_WT,highest_WT,lowest_NPA,highest_NPA)
a<-unique(top_FC)
plotheatmap(data.matrix(a))


