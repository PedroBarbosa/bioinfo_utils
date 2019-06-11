library(biomaRt)
library(DESeq2)
library("RColorBrewer")
library("apeglm")
library("pheatmap")
library(fgsea)
library(ggplot2)

explore_data_based_on_transformed_variance <- function(dds, cond){
  #Deseq2 offers 2 transformation methods to stabilize variance across the mean
  #Which method to choose, good discussion in the vignette
  #This is good explore data by  clustering samples by their distance of making a PCA analysis without the bias in the variance for high counts genes
  show("Stabilizing the variance to explore the data..")
  #v_transformed <- vst(dds, blind = FALSE)
  v_transformed <- rlog(dds, blind = FALSE) 
  #v_transformed <- log2(counts(dds) + 1)
  
  #PCA
  plot(DESeq2::plotPCA(v_transformed, intgroup = cond))

  #Sample distances
  sampleDists <- dist(t(assay(v_transformed)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(colData(v_transformed)[,cond])
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)

}


visualize_differences_in_variance_transformation <- function(log_t, vsd, rld) {
  
  #Show the effect of transformation plotting two samples against each other
  df <- bind_rows(
    as_data_frame(log_t[, 1:2]+1) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
  ggplot(df, aes(x = df$myC_time0_rep1, y = df$myC_time12_rep1)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 
  
}


##################################
######Differential testing########
##################################
#Recommended the DESeq function applied to raw counts,  which also takes into account the dependence of the variance
#of counts on the mean value during the dispersion estimation step.
run_DEseq_tests <- function(dds, test_contrast, coefficient, log2cutoff, padjcutoff) {
  dds_2 <- DESeq(dds)
  res <- results(dds_2, contrast=test_contrast, alpha=padjcutoff, test="Wald", independentFiltering = TRUE, pAdjustMethod="BH")
  summary(res)
  #res_24_0 <- results(dds, contrast=c("timepoints","24","0"), alpha=0.05, pAdjustMethod="BH")
  if (length(resultsNames(res)) > 0) {
    res <- lfcShrink(dds_2, coef = coefficient, type="apeglm", res= res)
  } else {
    res <- lfcShrink(dds_2, type="ashr", contrast = test_contrast)
  }

  res_sign <- subset(res, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
  return(list(res,res_sign))
  
}

sort_byFC <- function(res){
  sort_down <- res[ order(res$log2FoldChange), ]
  sort_up <- res[ order(res$log2FoldChange, decreasing = TRUE), ]
  return(list(sort_up,sort_down))
}


plot_single_gene_variation <- function(dds, group, gene){
  #Normalized data
  plotCounts(dds, gene = gene, intgroup=group)
  
  #From counts
  library("ggbeeswarm")
  geneCounts <- plotCounts(dds, gene = gene, intgroup = group,
                           returnData = TRUE)
  ggplot(geneCounts, aes(x = timepoints, y = count, color = group)) +
    scale_y_log10() +  geom_beeswarm(cex = 3)
}


plot_ma <- function(res_shrinked, padjcutoff){
  plotMA(res_shrinked, alpha=padjcutoff, xlab="mean of normalized counts", ylim = c(-5, 5))
  
  #with individual labels
  #topGene <- rownames(ma_12_0)[which.min(ma_12_0$padj)]
  #with(ma_12_0[topGene, ], {
  #  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  #  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
  #})
}

plot_volcano <- function(res_shrinked, padjcutoff, log2cutoff, title){
  res_shrinked$up <- res_shrinked$log2FoldChange > log2cutoff & res_shrinked$padj < padjcutoff
  res_shrinked$down <- res_shrinked$log2FoldChange < -log2cutoff & res_shrinked$padj < padjcutoff
  res_shrinked$threshold <- as.factor(abs(res_shrinked$log2FoldChange) > log2cutoff & res_shrinked$padj < padjcutoff)
  
  volcano_plot <- ggplot(data=res_shrinked, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(data=res_shrinked, size=1, colour="gray") +
    geom_point(data=res_shrinked[res_shrinked$down==TRUE, ], size=1, colour="#CC0000") +
    geom_point(data=res_shrinked[res_shrinked$up==TRUE, ], size=1, colour="#000099") +
    xlab("log2 fold change") +
    ylab("-log10 p-value adjusted") +
    ggtitle(title)+
    scale_y_continuous() +
    theme_bw() +
    #theme(axis.title.y = element_text(face="bold", size=16),
    #      axis.title.x = element_text(face="bold", size=16, colour="black"),
    #      axis.text = element_text(size=12),
    #      legend.title =element_blank() ,
    #      legend.text = element_text(size = 12)) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(volcano_plot)
}




cluster_sign_genes <- function(dds, sigGenes, group){
  print(paste("Number of DE genes across all combinations: ", length(sigGenes)))
  print("Clustering them based on the variance transformed raw counts.")
  v_transformed <- vst(dds, blind = FALSE)
  # print(typeof)
  #counts_justSigGenes <- counts(dds)[sigGenes,]
  counts_justSigGenes <- v_transformed[sigGenes,]
  DESeq2::plotPCA(counts_justSigGenes, intgroup = group)
}

annotate_results <- function(res, genome){
  if (genome == "mm10") {
    ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  }
  else if (genome == "hg38" | genome == "hg19" ) {
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  }
  
  genedesc <-  getBM(attributes=c('ensembl_gene_id_version','external_gene_name','description'), 
                     filters = 'ensembl_gene_id_version', values = rownames(res), mart =ensembl)
  
  res$symbol <- genedesc$external_gene_name[match(rownames(res),genedesc$ensembl_gene_id_version)]
  res$description <- genedesc$description[match(rownames(res),genedesc$ensembl_gene_id_version)]
  return(res)
}

#RuvSEQ specific
plot_umwanted_variation_factors <-function(pdata_normalized, groups, ksize) {
  par(mfrow = c(ksize, 1), mar = c(3,5,3,1))
  for (i in 1:ksize) {
    stripchart(pdata_normalized[, i] ~ groups, vertical = TRUE, main = paste0("W", i))
    abline(h = 0)
  }
  par(mfrow = c(1, 1))
}

run_fgsea_analysis <- function(de_table_annotated, hallmark_file, genesetname){
  #Note to GSEA users: Gene set enrichment analysis identifies gene sets consisting of co-regulated genes;
  #GO gene sets are based on ontologies and do not necessarily comprise co-regulated genes.
  de_table_annotated$stat <- de_table_annotated$log2FoldChange / de_table_annotated$lfcSE
  de_table_annotated <- as_tibble(de_table_annotated)
  res <- de_table_annotated %>% 
    dplyr::select(symbol, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(symbol) %>% 
    summarize(stat=mean(stat))
  
  
  ranks <- deframe(res)
  #pathways.reactome <- reactomePathways(names(ranks))
  pathways.hallmark <- gmtPathways(hallmark_file)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
  
  #plot(plotEnrichment(fgseaRes[["HALLMARK_MYOGENESIS"]], ranks))
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  #plot(plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, 
  #              gseaParam = 0.5))
  

  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    DT::datatable()
  
  if (nrow(fgseaResTidy) > 50){
    fgseaResTidy <- fgseaResTidy %>% top_n(50, wt=abs(-NES)) 
  }
  plot(ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Term", y="Normalized Enrichment Score",
         title=paste0(genesetname, " NES from GSEA")) + 
    theme_minimal())
  
  return(fgseaResTidy)
}

run_goseq_analysis <- function(de_table_annotated, genome, DE_genes){
  library(goseq)
  genes=as.integer(as.logical(rownames(de_table_annotated) %in% DE_genes))
  names(genes) = rownames(de_table_annotated)
  
  #genes=as.integer(de_table_annotated$padj[!is.na(de_table_annotated$padj)])
  #names(genes)=row.names(de_table_annotated[!is.na(de_table_annotated$padj),])
  names(genes) <- gsub("\\..*", "", names(genes))
  show(table(genes))
  supportedOrganisms()[supportedOrganisms()$Genome==genome,]
  pwf=nullp(genes,genome,"ensGene")
  
  go.results=goseq(pwf,genome,"ensGene", test.cats=c("GO:BP"), use_genes_without_cat=TRUE)
  plot(go.results %>% 
    top_n(10, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="p value", size="Count"))
  return(go.results)
}

run_analysis <- function(dds, group_combination, log2cutoff, padjcutoff, genome){
  library(ggplot2)
  #Data exploration
  explore_data_based_on_transformed_variance(dds, group_combination[[1]])
  
  
  #Diff expression
  out <-run_DEseq_tests(dds, group_combination,
                        paste(group_combination[[1]], group_combination[[2]], "_vs_", group_combination[[3]]), log2cutoff, padjcutoff)
  out_sign_annot <- annotate_results(out[[2]], genome)            
  outbasename <- paste(group_combination[[2]], "_vs_", group_combination[[3]], sep="")
  print(paste0("Number of genes in the ", outbasename, ": ", nrow(out[[2]])))
  write.table(out_sign_annot[,c("symbol","description","log2FoldChange","padj")], quote= FALSE, sep="\t", 
              file=paste0(outbasename,".csv", sep=""),row.names=TRUE)
  plot_ma(out[[1]] , padjcutoff)
  
  #Volcano
  plot_volcano(as.data.frame(out[[1]]), padjcutoff, log2cutoff, outbasename)
  
  #Plot one specific gene
  #topGene <- rownames(out[[1]])[which.min(out[[1]]$padj)]
  #plot_single_gene_variation(dds, group = "group_combination", topGene )
  return(list(out[[1]], out_sign_annot))
}

