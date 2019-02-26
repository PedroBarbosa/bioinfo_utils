library(ggplot2)
library(RColorBrewer)

setwd("C:/Users/PsBBa/Desktop")
df <- read.table("mappingStats_hg19_preprocessed_withPlateID.tsv", sep="\t",comment.char = "",check.names = FALSE,header=TRUE)

#ALIGNMENTS VS UNMAPPED
#principal component to color the direction of the correlations 
df$pc <- predict(prcomp(~df$`#alignments`+df$`#unmapped_reads`, df))[,1]

cols  <- brewer.pal(n=9, "Set1")
length((cols))
ggplot(df, aes(df$`#alignments`, df$`#unmapped`, color=df$plateID)) + 
 geom_point(shape=16,size=2) + 
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(axis.text=element_text(size=13),axis.title=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_y_continuous(name="Number of unmapped reads", labels = scales::comma) +
  scale_x_continuous(name="Number of aligments", label = scales::comma)


###BOXPLOTS###
fraction_proper=df$`#proper_pair`/df$`#alignments`*100
fraction_unique=df$`#unique`/df$`#alignments`*100
fraction_primary_q10=df$`#primary_linear_q10`/df$`#alignments`*100

fraction_unmapped=df$`#unmapped_reads`/df$`#reads`*100
fraction_duplicates=df$`#duplicate_reads`/df$`#reads`*100
fraction_chimeric=df$`#chimeric`/df$`#alignments`*100
fraction_seconday=df$`#secondary`/df$`#alignments`*100


labels_goodmetrics=c("primary_linear_q10","properly_mapped","uniquely_mapped")
labels_othermetrics=c("unmapped","duplicates","chimeric","secondary")

# get the greys (stolen from https://github.com/zonination/perceptions/blob/master/percept.R)
palette <- brewer.pal("Greys", n=9)
color.background = palette[2]
color.grid.major = palette[5]
cols  <- colorRampPalette(brewer.pal(8, "Set2"), alpha=TRUE)(ncol(df))

pdf("bamGoodMetrics.pdf",width=5, height = 3)
par(bty="o", bg=palette[1], mar=c(2,7,1.5,2))
boxplot(fraction_primary_q10, fraction_proper, fraction_unique, horizontal=TRUE,  lty=1,  boxwex=0.8,
        boxlwd=1, medlwd=1, outpch=19,col=cols,pars=list(outcol=cols,lwd=1,outbg=cols), yaxt="n", xaxt="n")

# plot gridlines
#for (i in seq(0,100,by=10)) {
#  lines(c(i,i), c(0,20), col=palette[4])
#}

#for (i in seq(1,17,by=1)) {
#  lines(c(-5,105), c(i,i), col=palette[4])
#}
axis(side=1, at=seq(0,100,by=5), col.axis=palette[7], cex.axis=0.8, lty=0 ,tick=NA, line=-1)
axis(side=2, at=1:3, col.axis=palette[7], cex.axis=0.8, lty=0, tick=NA, labels=labels_goodmetrics, las=1)
title(xlab="% of alignments", col.lab=palette[7], cex.lab=0.8, mgp=c(1,2,2))
dev.off()


pdf("bamBadMetrics.pdf", width=5, height =3)
par(bty="o", bg=palette[1], mar=c(2,7,1.5,2))
boxplot(fraction_unmapped, fraction_duplicates, fraction_chimeric,fraction_seconday, horizontal=TRUE,  lty=1,  boxwex=0.8,
        boxlwd=1, medlwd=1, outpch=19,col=cols,pars=list(outcol=cols,lwd=1,outbg=cols), yaxt="n", xaxt="n")
axis(side=1, at=seq(0,100,by=5), col.axis=palette[7], cex.axis=0.8, lty=0 ,tick=NA, line=-1)
axis(side=2, at=1:4, col.axis=palette[7], cex.axis=0.8, lty=0, tick=NA, labels=labels_othermetrics, las=1)
title(xlab="% of alignments/reads", col.lab=palette[7], cex.lab=0.8, mgp=c(1,2,2))
dev.off()

