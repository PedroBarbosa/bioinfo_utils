library(psichomics)

setwd("~/Downloads/GTE_v7/")
gtex <- loadGtexData(tissue = "Heart")
annot <- loadAnnotation(listSplicingAnnotations()[[1]])
psi <- quantifySplicing(annot, gtex$GTEx$`Junction quantification`, genes="MYBPC3")



#################################
#########Preparing the data #####
#################################
setwd("~/Downloads/christian/splicing/fromStar/new")
listSplicingAnnotations()[[3]]
fixInNamespace("prepareGeneQuantSTAR", pos="package:psichomics")
#setnames(table, colnames(table)[[2]], paste0("col", index))
prepareGeneQuant("GD1"="GD1ReadsPerGene.out.tab",
                 "GD2"="GD2ReadsPerGene.out.tab",
                 "GD3"="GD3ReadsPerGene.out.tab",
                 "S_3y_CBE"="S_3y_CBEReadsPerGene.out.tab",
                 "S_10y_CBE"="S_10y_CBEReadsPerGene.out.tab",
                 "S_11y_CBE"="S_11y_CBEReadsPerGene.out.tab",
                 "S_3y_CTRL"="S_3y_CTRLReadsPerGene.out.tab",
                 "S_10y_CTRL"="S_10y_CTRLReadsPerGene.out.tab",
                 "S_11_CTRL"="S_11y_CTRLReadsPerGene.out.tab",
                 strandedness = "stranded (reverse)")


prepareJunctionQuant("GD1SJ.out.tab",
                     "GD2SJ.out.tab",
                     "GD3SJ.out.tab",
                     "S_3y_CBESJ.out.tab",
                     "S_10y_CBESJ.out.tab",
                     "S_11y_CBESJ.out.tab",
                     "S_3y_CTRLSJ.out.tab",
                     "S_10y_CTRLSJ.out.tab",
                     "S_11y_CTRLSJ.out.tab"
                     )

#####################################
#########Loading the data ###########
#####################################
data <- loadLocalFiles("~/Downloads/christian/splicing/fromStar/new")
geneExpr <- data$Data$`Gene expression`
junctionQuant <- data$Data$`Junction quantification`

groups <- list("GD"=colnames(geneExpr)[1:3], 
               "CBE"=colnames(geneExpr)[4:6], 
                "CTRL"=colnames(geneExpr)[7:9])

groups <- list("3Y" = colnames(geneExpr)[c(4,7)],
               "10Y" = colnames(geneExpr)[c(5,8)])

human <- listSplicingAnnotations()[[3]]
annotation <- loadAnnotation(human)

####################################
########### Diff expression ########
####################################
geneExprFiltered <- geneExpr[rowSums(geneExpr)>10,]
geneExprNorm <- normaliseGeneExpression(geneExprFiltered, log2transform=TRUE)
pcaGE_all <- performPCA(t(geneExprNorm))
plotPCA(pcaGE_all, groups=groups)
genes <- rownames(geneExprNorm)

###################################
############ Splicing #############
###################################
psi <- quantifySplicing(annotation, junctionQuant, minReads=10)
events <- rownames(psi)
head(events,1)
head(psi)
pca <- performPCA(t(psi[,-c(1,2,3)]), missingValues = 3)
plotPCA(pca, groups = groups, loadings = FALSE)

plotDistribution(psi["AFE_17_-_48575043_48604480_48574014_HOXB6", ], groups)
table <- calculateLoadingsContribution(pca)
head(table,1)
(tmp  <- grep("ENSG00000164308", genes, value=TRUE))
(tmp  <- grep("SE_1_-_17915_17368_17233_17055_WASH7P", events, value=TRUE))
plotDistribution(geneExprNorm[tmp, ], groups, psi=FALSE)

cbe_ctrl <- groups[c("CBE", "CTRL")]
diffSplicing <- diffAnalyses(psi, cbe_ctrl, analyses = "wilcoxRankSum")
dim(diffSplicing)
deltaPSIthreshold <- abs(diffSplicing$`∆ Median`) > 0.5
delta_psi <- diffSplicing[which(deltaPSIthreshold), ]
rownames(delta_psi)

