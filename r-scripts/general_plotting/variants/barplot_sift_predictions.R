library(ggplot2)
library(RColorBrewer)
variants_in_dataset=1914
variants_withLargeRearrangements_Discarded=1812
distinctVariantLocations=1807
distinctVariantLocations_withSiftScores=893
distinctVariantLocations_withCaddScores=1602

v<-c(variants_in_dataset,variants_withLargeRearrangements_Discarded,distinctVariantLocations,
     distinctVariantLocations_withSiftScores,distinctVariantLocations_withCaddScores)

cols <- brewer.pal(n = 5, name = "PuBuGn")
names(v)<-c("dataset","rearrangements removed",
            "distinct variant locations","sift matches", "cadd13 matches")
pdf(file = "Dataset_manual_match.pdf", width = 6, height = 6,)
par(xpd = T, mar = par()$mar + c(0,0,0,7))
xx <- barplot(v, ylab= "Number of records",ylim=c(0,2010), col = cols,xaxt='n', ann=FALSE)
text(x = xx, y = v, label = v, pos = 3, cex = 0.8, col = "darkgreen")
legend(locator(1), legend = names(v),fill = cols)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

# Create data
data=data.frame(name=c("#dataset","#rearrangements removed","#distinct variant locations",
                       "#distinct locations with Sift Scores", "#distinct locations with cadd socres") , 
                value=v)
# Barplot
#ggplot(data, aes(x=name, y=value)) + geom_bar(stat = "identity", width = 0.4) +
#  theme(axis.title.x = element_blank()) + labs(y="Number of records") + ylim(0,2000)
                   