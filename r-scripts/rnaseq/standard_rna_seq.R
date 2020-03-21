library(biomaRt)
library(DESeq2)
library("RColorBrewer")
library("apeglm")
library("pheatmap")
library(fgsea)
library(ggplot2)
library(tidyverse)
library(goseq)

ensembl_genes <- read_tsv("/Users/pbarbosa/analysis/genome_utilities/hg38/mart_hg38.txt") 
names(ensembl_genes)[names(ensembl_genes) == "Gene name"] <- "symbol"
names(ensembl_genes)[names(ensembl_genes) == "Gene description"] <- "description"
names(ensembl_genes)[names(ensembl_genes) == "Gene stable ID"] <- "gene_id"
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
  if (length(resultsNames(res)) > 0) {
    #res <- lfcShrink(dds_2, type="apeglm", res= res)
    res <- lfcShrink(dds_2, coef = coefficient, type="apeglm", res= res)
  } else {
    #res <- lfcShrink(dds_2, coef = coefficient, type="apeglm", res= res)
    #res <- lfcShrink(dds_2, type="ashr", contrast = test_contrast)
    res <- lfcShrink(dds_2, type="normal", contrast = test_contrast)
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
  DESeq2::plotMA(res_shrinked, padjcutoff, xlab="mean of normalized counts",ylim = c(-5, 5))
  
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
  
  nup <- sum(res_shrinked$up, na.rm = TRUE)
  ndown <- sum(res_shrinked$down, na.rm = TRUE) 
  annotations <- data.frame(
    xpos = c(-2, 2),
    ypos =  c(9, 9),
    annotateText = c(paste0("N=",ndown),
                     paste0("N=", nup)),
    col = c("#CC0000", "#000099"),
    size = c(3,3),
    hjustvar = c(0, 0) ,
    vjustvar = c(0, 0))
  
  
  volcano_plot <- ggplot(data=res_shrinked, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(data=res_shrinked, size=1, colour="gray") +
    geom_point(data=res_shrinked[res_shrinked$down==TRUE, ], size=1, colour="#CC0000") +
    geom_point(data=res_shrinked[res_shrinked$up==TRUE, ], size=1, colour="#000099") +
    geom_text(data=annotations,aes(x=xpos,y=ypos, label=annotateText, colour=col, label.size = size,
                                   show.legend = FALSE)) +
    scale_colour_manual(values=c("#000099","#CC0000")) + 
    xlab("log2 fold change") +
    ylab("-log10 p-value adjusted") +
    ylim(0,7.5) +
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

run_fgsea_analysis <- function(de_table_annotated, hallmark_file, genesetname, is_deseq = T, is_psichomics = F, top_n_value = 50){
  #Note to GSEA users: Gene set enrichment analysis identifies gene sets consisting of co-regulated genes;
  #GO gene sets are based on ontologies and do not necessarily comprise co-regulated genes.
  if(is_deseq){
    de_table_annotated$stat <- de_table_annotated$log2FoldChange / de_table_annotated$lfcSE    
  }
  else if(is_psichomics){
    pvalueThreshold <- de_table_annotated$`p-value(BH_adjusted)` < 0.1
    de_table_annotated = de_table_annotated[pvalueThreshold, ]
    show(paste("Genes to be ranked by LogFC after setting a p-value adjusted cutoff:",nrow(de_table_annotated)))
    de_table_annotated = de_table_annotated %>% rename("symbol" = "Gene name", "stat" = "log2_Fold_change")
  }

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
  
  #plotEnrichment(fgseaRes[["KEGG_CELLCYLE"]], ranks)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, 
                gseaParam = 0.5)
  

  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    DT::datatable()
  
  if (nrow(fgseaResTidy) > top_n_value){
    fgseaResTidy <- fgseaResTidy %>% top_n(top_n_value, wt=abs(NES)) 
    #fgseaResTidy <- fgseaResTidy %>% top_n(top_n_value, wt=-padj) 
  }

 
  fgseaResTidy$sig = factor(fgseaResTidy$padj<0.05, levels = c(TRUE, FALSE))
  #fgseaResTidy <- fgseaResTidy_kegg %>% filter(sig == T)
  plot(ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj < 0.05)) +
    scale_fill_discrete(limits = c('FALSE', 'TRUE')) + 
    coord_flip() +
    labs(x="Term", y="Normalized Enrichment Score",
         title=paste0(genesetname, " NES from GSEA")) + 
    theme_minimal())
  
  return(fgseaResTidy)
}

run_goseq_analysis <- function(de_table_annotated, genome, DE_genes, is_deseq=T){

  if(is_deseq) de_table_annotated <- de_table_annotated %>% remove_rownames() %>% column_to_rownames(var='gene_id')
  genes=as.integer(as.logical(rownames(de_table_annotated) %in% DE_genes))
  names(genes) = rownames(de_table_annotated)   

  #genes=as.integer(de_table_annotated$padj[!is.na(de_table_annotated$padj)])
  #names(genes)=row.names(de_table_annotated[!is.na(de_table_annotated$padj),])
  names(genes) <- gsub("\\..*", "", names(genes))
  show(table(genes))
  supportedOrganisms()[supportedOrganisms()$Genome==genome,]
  pwf=nullp(genes,genome,"ensGene")
  
  go.results=goseq(pwf,genome,"ensGene", test.cats=c("GO:BP"), use_genes_without_cat=FALSE)
  plot(go.results %>% 
    top_n(30, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="Over represented GO term", colour="p value", size="Count"))
  
  plot(go.results %>% 
         top_n(30, wt=-under_represented_pvalue) %>%
         mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
         ggplot(aes(x=hitsPerc, 
                    y=term, 
                    colour=under_represented_pvalue, 
                    size=numDEInCat)) +
         geom_point() +
         expand_limits(x=0) +
         labs(x="Hits (%)", y="Under represented GO term", colour="p value", size="Count"))
  return(go.results)
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


run_analysis <- function(dds, group_combination, log2cutoff, padjcutoff, genome, explore_data=FALSE){
  #Data exploration
  if (explore_data == TRUE){
    show("Data exploration is set to TRUE..")
    explore_data_based_on_transformed_variance(dds, group_combination[[1]])  
  }

  #Diff expression
  show("Running DE tests and shrinking log2FC..")
  out <-run_DEseq_tests(dds, group_combination,
                        paste(group_combination[[1]], group_combination[[2]], "_vs_", group_combination[[3]]), log2cutoff, padjcutoff)

  #out_sign_annot <- annotate_results(out[[2]], genome)
  rownames(out[[2]]) <- gsub("\\..*", "", rownames(out[[2]]))
  out_sign_annot <- merge(as(out[[2]], "data.frame") , ensembl_genes, all.x = T, by.x=0, by.y="gene_id") 
  out_sign_annot <-out_sign_annot[!duplicated(out_sign_annot$Row.names),]
  outbasename <- paste(group_combination[[2]], "_vs_", group_combination[[3]], sep="")
  print(paste0("Number of genes in the ", outbasename, ": ", nrow(out[[2]])))
  write.table(out_sign_annot[,c("Row.names", "symbol","description","log2FoldChange","padj")], quote= FALSE, sep="\t", 
              file=paste0(outbasename,".csv", sep=""),row.names=FALSE)
  
  show("Generating MA Plot")
  plot_ma(out[[1]] , padjcutoff)
  
  show("Generating Volcano Plot!!")
  #Volcano
  plot_volcano(as.data.frame(out[[1]]), padjcutoff, log2cutoff, outbasename)
  
  #Plot one specific gene
  #topGene <- rownames(out[[1]])[which.min(out[[1]]$padj)]
  #plot_single_gene_variation(dds, group = "group_combination", topGene )
  return(list(out[[1]], out_sign_annot))
}

