library(psichomics)
#################################
#########Preparing the data #####
#################################
setwd("~/analysis/cajal_bodies/psichomics/")
listSplicingAnnotations()[[3]]

fixInNamespace("prepareGeneQuantSTAR", pos="package:psichomics")
#setnames(table, colnames(table)[[2]], paste0("col", index))

prepareGeneQuant("t0_rep1"="myC_time0_rep1_ReadsPerGene.out.tab",
                 "t0_rep2"="myC_time0_rep2_ReadsPerGene.out.tab",
                 "t0_rep3"="myC_time0_rep3_ReadsPerGene.out.tab",
                 "t12_rep1"="myC_time12_rep1_ReadsPerGene.out.tab",
                 "t12_rep2"="myC_time12_rep2_ReadsPerGene.out.tab",
                 "t12_rep3"="myC_time12_rep3_ReadsPerGene.out.tab",
                 "t24_rep1"="myC_time24_rep1_ReadsPerGene.out.tab",
                 "t24_rep2"="myC_time24_rep2_ReadsPerGene.out.tab",
                 "t24_rep3"="myC_time24_rep3_ReadsPerGene.out.tab",
                 "t48_rep1"="myC_time48_rep1_ReadsPerGene.out.tab",
                 "t48_rep2"="myC_time48_rep2_ReadsPerGene.out.tab",
                 "t48_rep3"="myC_time48_rep3_ReadsPerGene.out.tab",
                 strandedness = "stranded (reverse)")


prepareJunctionQuant("t0_rep1"="myC_time0_rep1_SJ.out.tab",
                     "t0_rep2"="myC_time0_rep2_SJ.out.tab",
                     "t0_rep3"="myC_time0_rep3_SJ.out.tab",
                     "t12_rep1"="myC_time12_rep1_SJ.out.tab",
                     "t12_rep2"="myC_time12_rep2_SJ.out.tab",
                     "t12_rep3"="myC_time12_rep3_SJ.out.tab",
                     "t24_rep1"="myC_time24_rep1_SJ.out.tab",
                     "t24_rep2"="myC_time24_rep2_SJ.out.tab",
                     "t24_rep3"="myC_time24_rep3_SJ.out.tab",
                     "t48_rep1"="myC_time48_rep1_SJ.out.tab",
                     "t48_rep2"="myC_time48_rep2_SJ.out.tab",
                     "t48_rep3"="myC_time48_rep3_SJ.out.tab"
                     )

#####################################
#########Loading the data ###########
#####################################
data <- loadLocalFiles(".")
geneExpr <- data$Data$`Gene expression`
junctionQuant <- data$Data$`Junction quantification`
colnames(geneExpr)
groups <- list("t0"=colnames(geneExpr)[1:3], 
               "t12"=colnames(geneExpr)[4:6], 
               "t24"=colnames(geneExpr)[7:9],
               "t48"=colnames(geneExpr)[10:12])

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

cbe_ctrl <- groups[c("t0", "t12")]
diffSplicing <- diffAnalyses(psi, cbe_ctrl, analyses = "wilcoxRankSum")
dim(diffSplicing)
deltaPSIthreshold <- abs(diffSplicing$`∆ Median`) > 0.5
delta_psi <- diffSplicing[which(deltaPSIthreshold), ]
rownames(delta_psi)
