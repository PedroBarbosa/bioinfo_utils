library(biomaRt)
library(DESeq2)
library("RColorBrewer")
library("apeglm")
library("pheatmap")
library(fgsea)
library(gprofiler2)
library(ggplot2)
library(tidyverse)
library(goseq)
library(writexl)

ensembl_genes <- read_tsv("/Users/pbarbosa/MEOCloud/analysis/genome_utilities/hg38/mart_hg38_v33.txt") 
names(ensembl_genes)[names(ensembl_genes) == "Gene name"] <- "symbol"
names(ensembl_genes)[names(ensembl_genes) == "Gene description"] <- "description"
names(ensembl_genes)[names(ensembl_genes) == "Gene stable ID"] <- "gene_id"

########################################################
################# DATA EXPLORATION #####################
########################################################
explore_data_based_on_transformed_variance <- function(dds, cond){
  #Deseq2 offers 2 transformation methods to stabilize variance across the mean
  #Which method to choose, good discussion in the vignette
  #This is good explore data by  clustering samples by their distance of making a PCA analysis without the bias in the variance for high counts genes
  show("Stabilizing the variance to explore the data..")
  #v_transformed <- vst(dds, blind = TRUE)
  v_transformed <- rlog(dds, blind = TRUE) 
  #v_transformed <- log2(counts(dds) + 1)
  
  #PCA
  plot(DESeq2::plotPCA(v_transformed, intgroup = cond))

  #Sample distances
  sampleDists <- dist(t(assay(v_transformed)))
  plot(hclust(dist(sampleDists)))
  sampleDistMatrix <- as.matrix( sampleDists )
  #rownames(sampleDistMatrix) <- paste(colData(v_transformed))#ß[,cond])
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


##################################################
############# DIFFERENTIAL TESTING ###############
##################################################
#Recommended the DESeq function applied to raw counts,  which also takes into account the dependence of the variance
#of counts on the mean value during the dispersion estimation step.

#Contrast vs coef 
#They are not identical. Using contrast is similar to what DESeq2 used to do: it forms an expanded model matrix,
#treating all factor levels equally, and averages over all distances between all pairs of factor levels to estimate the prior.
#Using coef, it just looks at that column of the model matrix (so usually that would be one level against the reference level)
#and estimates the prior for that coefficient from the distribution of those MLE of coefficients. 
#I implemented both for lfcShrink, because 'contrast' provides backward support (letting people get the same coefficient they obtained with previous versions), 
#while future types of shrinkage estimators will use the 'coef' approach, which is much simpler.
#However, for the 'coef' approach, the shrinkage depends on which level is chosen as reference.
#Expanded model matrices helped to avoid this, but introduced a lot of complexity into the modeling steps.
run_DEseq_tests <- function(dds, test_contrast, coefficient, log2cutoff, padjcutoff, use_contrast, shrinkage_method) {
  dds <- DESeq(dds)
  if (use_contrast){
    res <- results(dds, contrast = test_contrast, alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
    if (shrinkage_method == "apeglm") show("apeglm shrinkage method is not allowed when using contrast (instead of coeff)")
    res <- lfcShrink(dds, contrast = test_contrast, type=shrinkage_method, res = res)
  } else{
    show(resultsNames(dds))
    res <- results(dds, name = coefficient, alpha=padjcutoff, test="Wald", independentFiltering = T, pAdjustMethod="BH")
    res <- lfcShrink(dds, coef = coefficient, type=shrinkage_method, res= res)
  }
  summary(res)

  res_sign <- subset(res, padj<=padjcutoff & abs(log2FoldChange)>=log2cutoff)
  return(list(res,res_sign))
}

sort_byFC <- function(res){
  sort_down <- res[ order(res$log2FoldChange), ]
  sort_up <- res[ order(res$log2FoldChange, decreasing = TRUE), ]
  return(list(sort_up,sort_down))
}


###############################################
########### GENE EXPRESSION RESULTS ###########
###############################################
plot_single_gene_variation <- function(dds, group, gene){
  plotCounts(dds, gene = gene, intgroup=group)
}


plot_ma <- function(res_shrinked, padjcutoff){
  DESeq2::plotMA(res_shrinked, padjcutoff, xlab="Mean of normalized counts",ylim = c(-5, 5))
  #with individual labels
  #topGene <- rownames(ma_12_0)[which.min(ma_12_0$padj)]
  #with(ma_12_0[topGene, ], {
  #  points(baseMean, log2FoldChange, col="darkblue", cex=2, lwd=2)
  #  text(baseMean, log2FoldChange, topGene, pos=2, col="darkblue")
  #})
}

plot_volcano <- function(res_shrinked, padjcutoff, log2cutoff, title){
  show(paste0("Number of genes with less mean normalized counts than the optimal threshold (padj NA): ", nrow(res_shrinked[is.na(res_shrinked$padj),])))
  res_shrinked <- res_shrinked[!is.na(res_shrinked$padj),]
  res_shrinked$up <- res_shrinked$log2FoldChange > log2cutoff & res_shrinked$padj < padjcutoff
  res_shrinked$down <- res_shrinked$log2FoldChange < -log2cutoff & res_shrinked$padj < padjcutoff
  res_shrinked$threshold <- as.factor(abs(res_shrinked$log2FoldChange) > log2cutoff & res_shrinked$padj < padjcutoff)
  nup <- sum(res_shrinked$up, na.rm = TRUE)
  ndown <- sum(res_shrinked$down, na.rm = TRUE) 
  annotations <- data.frame(
    xpos = c(-3, 3),
    #ypos = c(max(-log10(res_shrinked$padj)) + 1, max(-log10(res_shrinked$padj))+ 1),
    ypos =  c(30,30),
    annotateText = c(paste0("N=",ndown),
                     paste0("N=",nup)),
    col = c("#CC0000", "#000099"),
    size = c(20,20),
    hjustvar = c(0, 0) ,
    vjustvar = c(0, 0))
  
  volcano_plot <- ggplot(data=res_shrinked, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(data=res_shrinked, size=1, colour="gray") +
    geom_point(data=res_shrinked[res_shrinked$down==TRUE, ], size=1, colour="#000099") +
    geom_point(data=res_shrinked[res_shrinked$up==TRUE, ], size=1, colour="#CC0000") +
    geom_text(data=annotations,aes(x=xpos,y=ypos, label=annotateText, colour=col)) +
    scale_colour_manual(values=c("#CC0000", "#000099")) + 
    theme_bw() +
    theme(panel.grid = element_blank(), text = element_text(size = 20)) +
    geom_hline(yintercept=-log10(padjcutoff), size = 0.3) + 
    geom_vline(xintercept = c(-log2cutoff, log2cutoff), size = 0.3) +
    xlab("Log Fold Change") +
    ylab("-Log10 p-value adjusted") +
    #coord_cartesian(xlim =c(min(res_shrinked$log2FoldChange) - 0.2, max(res_shrinked$log2FoldChange) + 0.2), 
    #                ylim =c(0, max(-log10(res_shrinked$padj) + 2)), expand = F) +
    #ylim(0, 17) +
    scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
    theme(legend.position="none")
    
  #theme(axis.title.y = element_text(face="bold", size=16),
    #      axis.title.x = element_text(face="bold", size=16, colour="black"),
    #      axis.text = element_text(size=12),
    #      legend.title =element_blank() ,
    #      legend.text = element_text(size = 12)) +

  plot(volcano_plot)
}

##############################
#### Expression heatmaps #####
##############################
plot_expression_heatmaps <- function(genes_expression_data, cluster_rows = T, cluster_cols = T, show_rownames = F, color_pallete = bluered(100)){
  #Using heatmap.2
  library("gplots")
  heatmap.2(genes_expression_data, Colv = cluster_cols, Rowv = cluster_rows, 
            labRow = show_rownames,
            col = color_pallete, 
            scale = "row",  
            trace = "none",
            density.info = "none", 
            cexCol = 1, 
            #keysize = 1,
            margins = c(5,8),
            srtCol=45)
  
  #Using pheatmap
  library(pheatmap)
  pheatmap(genes_expression_data, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
           show_rownames = show_rownames,
           col = color_pallete,
           scale = "row",
           angle_col = 45,
           fontsize = 25,
           legend = T
  )
}


#########################
#### GENE ANNOTATION ####
#########################
annotate_results <- function(res, annotate_locally, genome){
  res <- as.data.frame(res)
  rownames(res) <- gsub("\\..*", "", rownames(res))
  
  if (annotate_locally){
    out_sign_annot <- left_join(rownames_to_column(res), 
                                dplyr::select(ensembl_genes, c(gene_id, symbol, description)), 
                                by = c( "rowname" = "gene_id"))
    out_sign_annot <-out_sign_annot[!duplicated(out_sign_annot$rowname),] %>% remove_rownames %>% column_to_rownames()
    return(out_sign_annot)
  }
  
  
  else{
    if (genome == "mm10") {
      ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
    }
    else if (genome == "hg38" | genome == "hg19" ) {
      ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    }
    genedesc <-  getBM(attributes=c('ensembl_gene_id','external_gene_name','description'), 
                       filters = 'ensembl_gene_id', values = rownames(res), mart =ensembl)
    res$symbol <- genedesc$external_gene_name[match(rownames(res),genedesc$ensembl_gene_id)]
    res$description <- genedesc$description[match(rownames(res),genedesc$ensembl_gene_id)]
    return(res)
  }
}



################################################
########## ENRICHMENT ANALYSIS #################
################################################
run_enrichment_gprofiler <- function(gene_list, custom_genes_background=NULL, organism=c("hsapiens", "mmusculus"),
                                     ordered_gene_query = F, multi_query = F, retrieve_only_sign_results = T, measure_under = F,
                                     evcodes = F, correction_method = c("gSCS", "fdr", "bonferroni"),
                                     exclude_iea = T, domain_scope = c("annotated", "known"), sources = NULL, retrieve_short_link = F){
  organism <- match.arg(organism)
  correction_method <- match.arg(correction_method)
  domain_scope <- match.arg(domain_scope)
  if (!is.null(custom_genes_background)) { 
    domain_scope = "custom_annotated"
  }
  if (!is.null(sources)) { 
    sources = sources
  }
  else{
    sources = c("GO", "KEGG", "REAC")
  }
 

  gostres <- gost(query = gene_list,
                  organism = organism,
                  ordered_query = ordered_gene_query,
                  multi_query = multi_query,
                  significant = retrieve_only_sign_results,
                  exclude_iea = exclude_iea,
                  measure_underrepresentation = measure_under,
                  evcodes = F,
                  user_threshold = 0.05,
                  correction_method = correction_method,
                  domain_scope = domain_scope,
                  custom_bg = custom_genes_background,
                  numeric_ns = "",
                  sources = sources,
                  as_short_link = retrieve_short_link)
  
  if (is.null(gostres)){
    stop(paste0("Results are empty. It appears there is no significant terms within the sources tested: ", paste(sources,collapse=";")))
  }
  else if (retrieve_short_link){
    return(gostres)
  }
  
  show(gostplot(gostres))
  if (multi_query){
    cols_to_extract <- c("term_id", "term_name", "p_values", "term_size", "query_sizes", "intersection_sizes", "effective_domain_size", "source", "significant")
    results_table <- gostres$result[,cols_to_extract]
  } else{
    cols_to_extract <- c("term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size", "effective_domain_size", "source", "significant")
    results_table <- gostres$result[,cols_to_extract]
    results_table$ratio_present <- results_table$intersection_size / results_table$term_size
    results_table$log10padj <- -log10(results_table$p_value)
  }
  return(results_table)
}


fgsea_function <- function(pathways, ranks, top_n, npermutations, source, outlist){

  fgseaRes <- fgsea(pathways=pathways, 
                    stats=ranks, 
                    nperm=npermutations,
                    minSize = 10)
  
  #plotEnrichment(fgseaRes[["KEGG_CELLCYLE"]], ranks)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=7), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=7), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
                gseaParam = 0.5, render = F)
  
  fgseaRes$source = source
  fgseaRes$log10padj <- -log10(fgseaRes$padj)
  fgseaRes$significant = factor(fgseaRes$padj < 0.05, levels = c(TRUE, FALSE))
  fgseaRes <- dplyr::rename(fgseaRes, c(term_name = pathway, term_size = size)) 
  #term_size refers to the size of the pathway after removing genes from the ranked list
  #not present in the names of the pathway
  
  plot(top_n(fgseaRes, top_n, wt=log10padj) %>%
    top_n(top_n, wt=abs(NES)) %>%
    ggplot(aes(reorder(term_name, NES), NES)) +
         geom_col(aes(fill = padj < 0.05)) +
         scale_fill_discrete(limits = c('FALSE', 'TRUE')) + 
         coord_flip() +
         labs(x="Term", y="Normalized Enrichment Score",
              title=paste0(source, " NES from GSEA")) + 
         theme_minimal())
  
  outlist[[source]] = fgseaRes %>% dplyr::select(-leadingEdge)
  return(outlist)
}


run_fgsea_analysis <- function(gene_list, sources = NULL, is_ranked_already = F, 
                               gene_list_from = c("DESeq2", "psichomics", "rmats", "vast-tools", "majiq"), 
                               top_n = 25, npermutations = 1000){
  #Note to GSEA users: Gene set enrichment analysis identifies gene sets consisting of co-regulated genes;
  #GO gene sets are based on ontologies and do not necessarily comprise co-regulated genes.
  HALLMARKS_GMT <- "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/GSEA/h.all.v7.0.symbols.gmt"
  REAC_GMT <- "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/GSEA/c2.cp.reactome.v7.0.symbols.gmt"
  KEGG_GMT <- "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/GSEA/c2.cp.kegg.v7.0.symbols.gmt"
  GO_BP_GMT <- "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/GSEA/c5.bp.v7.0.symbols.gmt"
  GO_MF_GMT <- "/Users/pbarbosa/MEOCloud/analysis/genome_utilities/GSEA/c5.mf.v7.0.symbols.gmt"
  

  gene_list_from <- match.arg(gene_list_from)
  if (!is.null(sources)) { 
    sources <- sources
  }
  else{
    sources <- c("Hallmarks", "KEGG", "REAC", "GO:BP", "GO:MF")
  }
  
  if (!is_ranked_already){
    if(gene_list_from == "psichomics"){
      gene_list$stat <- gene_list$log2_Fold_change / gene_list$`p-value(BH_adjusted)`
      gene_list = gene_list %>% rename("symbol" = "Gene name")
    }
    else if(gene_list_from == "DESeq2"){
      gene_list$stat = gene_list$log2FoldChange / gene_list$lfcSE
    }
    
    else if (gene_list_from == "rmats" | gene_list_from == "vast-tools" | gene_list_from == "majiq"){
      gene_list <- gene_list  %>% mutate(stat = rank)
    }
    
    res <- as_tibble(gene_list) %>% 
      dplyr::select(symbol, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      dplyr::summarize(stat=mean(stat))
    ranks <- deframe(res)

  } else{
    ranks <- gene_list
  }
  
  outlist = list()
  for (source in sources){
    if (source == "Hallmarks"){
      pathways <- gmtPathways(HALLMARKS_GMT)
    }
    else if (source == "KEGG"){
      pathways <- gmtPathways(KEGG_GMT)
    }
    else if (source == "REAC"){
      pathways <- gmtPathways(REAC_GMT)
    }
    else if (source == "GO:BP"){
      pathways <- gmtPathways(GO_BP_GMT)
    }
    else if (source == "GO:MF"){
      pathways <- gmtPathways(GO_MF_GMT)
    }
    else{
      stop(paste0(source, " is not a valid source."))
    }
    outlist <- fgsea_function(pathways, ranks, top_n, npermutations, source, outlist)
  }
  out_df <- do.call("rbind", outlist)
  return(out_df)
}


run_goseq_analysis <- function(de_genes, all_genes, genome, sources=NULL, use_genes_without_cat = F, top_n = 15){
  
  if (!is.null(sources)) { 
    sources <- sources
  }
  else{
    sources <- c("GO:BP", "GO:MF", "GO:CC")
  }
  
  genes=as.integer(as.logical(all_genes %in% de_genes))
  names(genes) = all_genes 
  #genes=as.integer(de_table_annotated$padj[!is.na(de_table_annotated$padj)])
  #names(genes)=row.names(de_table_annotated[!is.na(de_table_annotated$padj),])

  show(table(genes))
  supportedOrganisms()[supportedOrganisms()$Genome==genome,]
  pwf=nullp(genes,genome,"ensGene")
  
  outlist = list()
  for (source in sources) { 
    show(paste("Analysising GO category:", source))
    if (!source %in% c("GO:BP", "GO:MF", "GO:CC")) stop(paste0("Source is not valid: ", source))
    
    goseq.results  <- goseq(pwf, genome, id = "ensGene", test.cats = source, use_genes_without_cat = use_genes_without_cat)
    goseq.results <- dplyr::rename(goseq.results, c(term_name = term, source= ontology, term_size = numInCat, 
                                                    intersection_size = numDEInCat, term_id = category, 
                                                    pvalue_over = over_represented_pvalue, pvalue_under = under_represented_pvalue))
    
    goseq.results$pvalue_adj_over <- p.adjust(goseq.results$pvalue_over, method = "fdr")
    goseq.results$pvalue_adj_under <- p.adjust(goseq.results$pvalue_under, method = "fdr")
    goseq.results$ratio_present <- goseq.results$intersection_size / goseq.results$term_size
    goseq.results$significant <- factor(goseq.results$pvalue_adj_over < 0.05 | goseq.results$pvalue_adj_under < 0.05, levels = c(TRUE, FALSE))
    goseq.results$log10padj_over <- -log10(goseq.results$pvalue_adj_over)
    goseq.results$log10padj_under <- -log10(goseq.results$pvalue_adj_under)
    plot_goseq_results(goseq.results, top_n, goseq.results$source)
    outlist[[source]] = goseq.results
  }
  out_df <- do.call("rbind", outlist)
  return(out_df)
}


plot_goseq_results <- function(go.results, top_n, ontology){
  #Over
  results_over <- head(go.results[order(go.results$pvalue_adj_over),], top_n)
  results_under <- head(go.results[order(go.results$pvalue_adj_under),], top_n)
  
  plot(results_over %>% 
         #top_n(top_n, wt=-pvalue_adj_over) %>% 
         mutate(hitsPerc=intersection_size*100/term_size) %>% 
         ggplot(aes(x=hitsPerc, 
                    y=term_name, 
                    colour=pvalue_adj_over, 
                    size=intersection_size)) +
         geom_point() +
         theme( text = element_text(size = 15)) + 
         expand_limits(x=0) +
         labs(x="Hits (%)", y=paste0("Over represented GOs in ", ontology, " ontology"), colour="padj", size="Count")) 
        
  
  #Under
  plot(results_under %>% 
         #top_n(top_n, wt=-pvalue_adj_under) %>% 
         mutate(hitsPerc=intersection_size*100/term_size) %>% 
         ggplot(aes(x=hitsPerc, 
                    y=term_name, 
                    colour=pvalue_adj_under, 
                    size=intersection_size)) +
         geom_point() +
         theme( text = element_text(size = 15)) +
         expand_limits(x=0) +
         labs(x="Hits (%)", y=paste0("Under represented GOs in ", ontology, " ontology"), colour="padj", size="Count")) 
}


#gProfiler
plot_enrichment_results <- function(results_table, top_n = 15, rank_by = "log10padj", label = "-log10(corrected P-value)", 
                                    short_term_size = F, size_of_short_term = 100, add_info_to_labels = T, font_size = 20,
                                    reverse_order = F){
  if (short_term_size) results_table <- dplyr::filter(results_table, term_size < size_of_short_term)
  
  if (reverse_order){
    results_table <- head(results_table[order(results_table[rank_by]),], top_n)
  } else{
    results_table <- head(results_table[order(-results_table[rank_by]),], top_n)
  }

  ggplot(results_table, aes(x=reorder(term_name, get(rank_by)), y=get(rank_by))) +
    geom_bar(stat="identity", width = 0.4, color = "black") +
    labs(x = ""  , y = label) +
    scale_y_continuous(position = "right") +
    coord_flip() +
    {if(add_info_to_labels)scale_x_discrete(labels = paste0(rev(results_table$term_name), " (", rev(results_table$source), "), N=", rev(results_table$term_size)))} +
    theme_light() +
    #theme(aspect.ratio = 1/1.5) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = font_size))
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


#############################################
########## MAIN RUNNING FUNCTION ############
#############################################
run_analysis <- function(dds, group_combination, log2cutoff, padjcutoff, use_contrast = F,
                         FC_shrinkage_method = c("apeglm", "ashr", "normal"), 
                         genome = "hg38",
                         explore_data= F,
                         annotate_locally = F){
  
  if (explore_data == TRUE){
    show("Data exploration is set to TRUE..")
    explore_data_based_on_transformed_variance(dds, group_combination[[1]])  
  }
  else {
    FC_shrinkage_method <- match.arg(FC_shrinkage_method)
    outbasename <- paste(group_combination[[2]], "_vs_", group_combination[[3]], sep="")
    show("Running DE tests and shrinking log2FC..")
    out <-run_DEseq_tests(dds, group_combination,
                          paste0(group_combination[[1]], "_", group_combination[[2]], "_vs_", group_combination[[3]]), 
                          log2cutoff, 
                          padjcutoff,
                          use_contrast,
                          FC_shrinkage_method)

    show("Annotating DE genes..")
    out_sign_annot <- annotate_results(out[[2]], annotate_locally, genome)
    out_sign_annot %>% rownames_to_column(var = "gene_id") %>% 
      dplyr::select(c(gene_id, symbol, description, log2FoldChange, pvalue, padj)) %>%
      arrange(-log2FoldChange) %>%
      write_xlsx(., path = paste0(outbasename,".xlsx", sep=""),col_names=T, format_headers = T)
    
    show(paste0("Number of genes in the ", outbasename, ": ", nrow(out_sign_annot)))
    show("Generating MA Plot")
    plot_ma(out[[1]] , padjcutoff)
    
    show("Generating Volcano Plot!!")
    plot_volcano(as.data.frame(out[[1]]), padjcutoff, log2cutoff, outbasename)
    
    #Plot one specific gene
    #topGene <- rownames(out[[1]])[which.min(out[[1]]$padj)]
    #plot_single_gene_variation(dds, group = "groups", topGene )
    return(list(out[[1]], out_sign_annot))
  }
}

