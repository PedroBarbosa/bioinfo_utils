library("argparse")
parser <- ArgumentParser(description='Process gene lists/ranks to do functional enrichment of differential splicing analysis')

parser$add_argument('tool', nargs='+', choices=c("majiq", "vastools", "rmats"), help='Tools to analyse. If more than 1 tool is given \\
                    a combined approach (genes from each tool added) will also be performed (in addition to individual enrichment test per tool)')
parser$add_argument('-o', '--outbasename', required=T,help='Basename to write output files')
parser$add_argument('-s', '--species', choices=c("hsapiens", "mmusculus"), default="hsapiens", help='Species name')
parser$add_argument('-n', '--no_custom', action="store_true", help='Do not use custom set of background genes in the statistical test. \\
                    Default: Use custom set of genes if "--*_negative_genes" is set. In the combined approach, custom set of background \\
                    genes will be used only if all the tools supplied have its "--*negative_genes" set. If not, custom set wont be used, \\
                    even if this argument is not set')
parser$add_argument('-p', '--plot_short_term', action="store_true", help='Whether to plot a subset of significant terms that have a short size')
parser$add_argument('-t', '--term_size', type="integer", default=400, help='Size of a term to be considered as short by "--plot_short_term"')

parser$add_argument('--vastools_positive_genes', help='Positive genes for GO enrichment produced by vastools')
parser$add_argument('--vastools_negative_genes', default=NULL, help='vastools background genes for GO enrichment. \\
                    If not set, gprofiler will run without custom backround list.')
parser$add_argument('--vastools_ranks', help='vastools ranks for GSEA analysis')
parser$add_argument('--vastools_ranks_groups', help='Partial strings representing 2 groups (split by ","). These groups \\
            must match colnames in the tidy table passed in the "--vastools_ranks" arg. e.g. (mock, CBE) and will be used \\
            to calculate dPSI differences and thus rank genes for GSEA analysis')
parser$add_argument('--rmats_positive_genes', help='Positive genes for GO enrichment produced by rMATS')
parser$add_argument('--rmats_negative_genes', default=NULL, help='rMATS background genes for GO enrichment. \\
                    If not set, gprofiler will run without custom backround list.')
parser$add_argument('--rmats_ranks', help='rMATs ranks for GSEA analysis')

parser$add_argument('--majiq_positive_genes', help='Positive genes for GO enrichment produced by MAJIQ')
parser$add_argument('--majiq_negative_genes', default=NULL, help='MAJIQ background genes for GO enrichment. \\
                    If not set, gprofiler will run without custom backround list.')

args <- parser$parse_args()
USE_SHORT_TERMS <<- args$plot_short_term
SHORT_TERM_SIZE <<- args$term_size
ORGANISM <<- args$species
OUTBASENAME <<- args$outbasename

repos <- "http://cran.us.r-project.org"
packages <- list(
  list(name="BiocManager"),
  list(name="dplyr"),
  list(name="tidyr"),
  list(name="fgsea"),
  list(name="readr"),
  list(name="gprofiler2"),
  list(name="ggplot2"),
  list(name="biomaRt"),
  list(name="msigdbr")
)

for (i in 1:length(packages)) {
  package <- packages[[i]]
  installed <- require(package$name, character.only=TRUE)
  if (installed) {
    next
  }
  if (package$name == "fgsea"){
    BiocManager::install("fgsea") 
  } else {
    install.packages(package$name, repos=repos)
  }
  installed <- require(package$name, character.only=TRUE)
  if (!installed) {
    stop(paste("could not install", package$name))
  }
}


#' @title GO enrichment of splicing events
#' @description Given a character vector of positive genes (with any relevant splicing event) 
#' performs gene list functional enrichment analysis using gprofiler or goseq packages.
#' 
#' If a set of background set of genes is provided, it will serve as a custom set of genes to 
#' fetch terms that will compose the statistical domain to perform the enrichment tests.
#' 
#' If a list of character vectors is provided in the \code{gene_list} argument, a comparison
#' between gene lists will be performed by passing the \code{multi_query} argument to gprofiler.
#' It's useful to compare functional enrichment profiles between splicing events detected by
#' different tools.
#'  
#' The background list represents genes that meet similar read coverage criteria as the positive
#' gene list, so GO enrichment of highly expressed genes is avoided.
#'
#' @param gene_list a character vector (or a list of character vectors) containing gene IDs that contain 
#' differentially spliced events.
#' @param custom_genes_background a character vector containing gene IDs to serve as the background in the statistical
#' procedure (generated using vast-tools compare with --GO option set.). If not set, gprofiler will be used in 
#' the \code{package} argument because it can automatically fetch annotated lists of IDs for any given data
#' sources to enrich. As opposed to \code{positive_genes} this argument doesn't accept a list of several backrgound 
#' gene sets. Therefore, be careful on the choice of the background gene set when providing multiple positive gene lists.
#' @param ... Arbitary number of additional arguments to be passed to the enrichment functions (check gprofiler and goseq 
#' enrichment functions to know which arguments can be passed). \code{run_enrichment_gprofiler}
#' 
#' @return a dataframe with the GO enriched terms
run_enrichment_gprofiler <- function(gene_list, custom_genes_background=NULL, 
                                     organism=c("hsapiens", "mmusculus"),
                                     ordered_gene_query = F,
                                     multi_query = F, 
                                     retrieve_only_sign_results = T, 
                                     measure_under = F,
                                     evcodes = F,
                                     correction_method = c("gSCS", "fdr", "bonferroni"),
                                     exclude_iea = T, 
                                     domain_scope = c("annotated", "known"), 
                                     sources = NULL,
                                     retrieve_short_link = F){
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
    show(paste0("Results are empty. It appears there is no significant terms within the sources tested: ", paste(sources,collapse=";")))
    return(NULL)
  }
  else if (retrieve_short_link){
    return(gostres)
  }
  
  if (multi_query){
    cols_to_extract <- c("term_id", "term_name", "p_values", "term_size", "query_sizes", "intersection_sizes", "effective_domain_size", "source", "significant")
    results_table <- gostres$result[,cols_to_extract]
  } else{
    cols_to_extract <- c("term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size", "effective_domain_size", "source", "significant")
    results_table <- gostres$result[,cols_to_extract]
    results_table$ratio_present <- results_table$intersection_size / results_table$term_size
    results_table$log10padj <- -log10(results_table$p_value)
  }
  return(list(results_table, gostres))
}


fgsea_function <- function(pathways, ranks, top_n, npermutations, source, outlist){
  
  fgseaRes <- fgseaMultilevel(pathways=pathways, 
                    stats=ranks, 
                    eps=0,
                    minSize = 5)
  
  #plotEnrichment(fgseaRes[["KEGG_CELLCYLE"]], ranks)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=7), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=7), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
                gseaParam = 0.5, render = T)
  
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
                               gene_list_from = c("psichomics", "rmats", "vast-tools", "majiq"), 
                               top_n = 25, npermutations = 2000, organism=c("hsapiens","mmusculus")){

  gene_list_from <- match.arg(gene_list_from)
  organism <- match.arg(organism)
  
  if (organism == "mmusculus"){
    show("Retrieving all mouse gene sets in MSigDB")
    m_df = msigdbr(species = "Mus musculus")
  } else if (organism == "hsapiens"){
    show("Retrieving all human gene sets in MSigDB")
    m_df = msigdbr(species = "Homo sapiens")
  }

  Hallmarks <- m_df %>% dplyr::filter(gs_cat == "H") %>% split(x = .$gene_symbol, f = .$gs_name)
  Reactome <- m_df %>% dplyr::filter(gs_subcat == "CP:REACTOME") %>% split(x = .$gene_symbol, f = .$gs_name)
  Kegg <- m_df %>% dplyr::filter(gs_subcat == "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
  GO_BP <- m_df %>% dplyr::filter(gs_subcat == "BP") %>% split(x = .$gene_symbol, f = .$gs_name)
  GO_MF <-  m_df %>% dplyr::filter(gs_subcat == "MF") %>% split(x = .$gene_symbol, f = .$gs_name)
  #gene_list <- gene_list %>% rename(symbol_mouse = symbol)
  #gene_list %>% left_join(dplyr::select(human_gene_symbol, gene_symbol) %>% distinct(), by=c("symbol_mouse"="gene_symbol")) %>%
  #  rename(symbol = human_gene_symbol)
  
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

    else if (gene_list_from == "rmats" | gene_list_from == "vast-tools" | gene_list_from == "majiq"){
      gene_list <- gene_list  %>% mutate(stat = rank)
    }

    res <- as_tibble(gene_list) %>% 
      dplyr::select(gene_id, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(gene_id) %>% 
      dplyr::summarize(stat=mean(stat))
   
    ranks <- tibble::deframe(res)
    
  } else{
    ranks <- gene_list
  }
  
  outlist = list()
  for (source in sources){
    if (source == "Hallmarks"){
      pathways <- Hallmarks
    }
    else if (source == "KEGG"){
      pathways <- Kegg
    }
    else if (source == "REAC"){
      pathways <- Reactome
    }
    else if (source == "GO:BP"){
      pathways <- GO_BP
    }
    else if (source == "GO:MF"){
      pathways <- GO_MF
    }
    else{
      stop(paste0(source, " is not a valid source."))
    }
    outlist <- fgsea_function(pathways, ranks, top_n, npermutations, source, outlist)
  }
  out_df <- do.call("rbind", outlist)
  return(out_df)
}


plot_enrichment_results <- function(results_table, top_n = 15, rank_by = "log10padj", label = "-log10(corrected P-value)", 
                                    short_term_size = F, size_of_short_term = 400, add_info_to_labels = T, font_size = 12,
                                    reverse_order = F){
  results_table <- as_tibble(results_table)
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

process_vastools <- function(pos, neg, ranks, ranks_groups, no_custom){
  #####################
  ### GO enrichment ###
  #####################
  if (is.null(pos)){
    stop("vastools positive genes file is missing. Please set '--vastools_positive_genes'")
  }
  vast_pos <- read_tsv(pos, col_names = "gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
  
  if (is.null(neg) || no_custom){
    background_genes <- NULL
  } else {
    background_genes <- read_tsv(neg, col_names="gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
  }
  
  vast_res <- run_enrichment_gprofiler(vast_pos, custom_genes_background = background_genes,
                                     retrieve_only_sign_results = T,
                                     exclude_iea = F,
                                     organism = ORGANISM,
                                     retrieve_short_link = F,
                                     measure_under = F, 
                                     domain_scope = "annotated")

  if(!is.null(vast_res)){
    pdf(paste0(OUTBASENAME, "_gostplot_vastools.pdf"))
    p <- gostplot(vast_res[[2]], interactive = FALSE)
    plot(p)
    dev.off()
    
    pdf(paste0(OUTBASENAME, "_gprofiler_vastools.pdf"))
    p <- plot_enrichment_results(vast_res[[1]], short_term_size=USE_SHORT_TERMS, size_of_short_term=SHORT_TERM_SIZE)
    plot(p)
    dev.off()
    write.table(vast_res[[1]], file = paste0(OUTBASENAME, "_gprofiler_vastools.tsv"), row.names = F, quote=F, sep="\t")
  }
  
  
  #####################
  ###### FGSEA ########
  #####################
  if (!is.null(ranks)){
    if(is.null(ranks_groups)){
      stop("If '--vastools_ranks' is set, you must pass '--vastools_ranks_groups' as well.")
    }
    
    to_rank <- read_tsv(ranks) %>% tidyr::separate(col = EVENT, into = c("symbol", "eventid"), sep="=")
    groups <- strsplit(ranks_groups, split=",")
    if (length(groups != 2)){
      stop("Two groups must exist in the '--vastools_ranks_groups' split by ','")
    }

    to_rank$groups[[1]] <- rowMeans(subset(to_rank, select = grepl(groups[[1]], colnames(to_rank))), na.rm = T)
    to_rank$groups[[2]] <- rowMeans(subset(to_rank, select = grepl(groups[[2]], colnames(to_rank))), na.rm = T)
    show(to_rank)
    
    ranks_vasttools <- to_rank %>% mutate(dPSI = groups[[2]] - groups[[1]]) %>% group_by(symbol) %>% 
      drop_na(dPSI) %>%
      summarise(rank=max(dPSI)) %>%
      filter(symbol != "") %>%
      distinct() %>%
      arrange(-rank)
    
    vast_res_gsea <- run_fgsea_analysis(ranks_vasttools, is_ranked_already = F, gene_list_from = "vast-tools", organism=organism)
  }
  return(list(vast_pos, background_genes))
}

process_majiq <- function(pos, neg, no_custom){
  #####################
  ### GO enrichment ###
  #####################
  if (is.null(pos)){
    stop("majiq positive genes file is missing. Please set '--majiq_positive_genes'")
  }
  majiq_pos <- read_tsv(pos, col_names = "gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
  
  if (is.null(neg) || no_custom){
    background_genes <- NULL
  } else {
    background_genes <- read_tsv(neg, col_names="gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
  }
  
  majiq_res <- run_enrichment_gprofiler(majiq_pos, custom_genes_background = background_genes,
                                       retrieve_only_sign_results = T,
                                       exclude_iea = F,
                                       organism = ORGANISM,
                                       retrieve_short_link = F,
                                       measure_under = F, 
                                       domain_scope = "annotated")
  
  if (!is.null(majiq_res)){
    pdf(paste0(OUTBASENAME, "_gostplot_majiq.pdf"))
    p <- gostplot(majiq_res[[2]], interactive=F)
    plot(p)
    dev.off()

    pdf(paste0(OUTBASENAME, "_gprofiler_majiq.pdf"))
    p <- plot_enrichment_results(majiq_res[[1]], short_term_size=USE_SHORT_TERMS, size_of_short_term=SHORT_TERM_SIZE)
    plot(p)
    dev.off()
    write.table(majiq_res[[1]], file = paste0(OUTBASENAME, "_gprofiler_majiq.tsv"), row.names = F, quote=F, sep="\t")
  } 
  return(list(majiq_pos, background_genes))
}
 
 
process_rmats <- function(pos, neg, ranks, no_custom){
    #####################
    ### GO enrichment ###
    #####################
    if (is.null(pos)){
      stop("rmats positive genes file is missing. Please set '--rmats_positive_genes'")
    }
    rmats_pos <- read_tsv(pos, col_names = "gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
    rmats_pos <- gsub("\\..*","", rmats_pos)
    
    if (is.null(neg) || no_custom){
      background_genes <- NULL
    } else {
      background_genes <- read_tsv(neg, col_names="gene_id") %>% filter(grepl("ENS", gene_id)) %>% pull()
      background_genes <- gsub("\\..*","", background_genes)
    }
    
    rmats_res <- run_enrichment_gprofiler(rmats_pos, custom_genes_background = background_genes,
                                         retrieve_only_sign_results = T,
                                         exclude_iea = F,
                                         organism = ORGANISM,
                                         retrieve_short_link = F,
                                         measure_under = F, 
                                         domain_scope = "annotated")
    if (!is.null(rmats_res)){
      pdf(paste0(OUTBASENAME, "_gostplot_rmats.pdf"))
      p <- gostplot(rmats_res[[2]], interactive = F)
      plot(p)
      dev.off()
      
      pdf(paste0(OUTBASENAME, "_gprofiler_rmats.pdf"))
      p <- plot_enrichment_results(rmats_res[[1]], short_term_size=USE_SHORT_TERMS, size_of_short_term=SHORT_TERM_SIZE)
      plot(p)
      dev.off()
      write.table(rmats_res[[1]], file = paste0(OUTBASENAME, "_gprofiler_rmats.tsv"), row.names = F, quote=F, sep="\t")
    }
    #####################
    ###### FGSEA ########
    #####################
    if (!is.null(ranks)){
      if (ORGANISM == "mmusculus") {
        ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
      }
      else{
        ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      }
  
      to_rank <- read_tsv(ranks, col_names = c("gene_id", "padj", "dPSI"))
      to_rank$gene_id <- gsub("\\..*","", to_rank$gene_id)
      
      genedesc <-  getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
                         filters = 'ensembl_gene_id', values = to_rank$gene_id , mart =ensembl)
      to_rank$symbol <- genedesc$external_gene_name[match(to_rank$gene_id, genedesc$ensembl_gene_id)]
      
      
      ranks_rmats <- to_rank %>% mutate(padj = ifelse(padj == 1, 0.9999999, ifelse(padj == 0, 1e-20, padj)),
                                        dPSI = -dPSI, value = -log10(padj) * dPSI) %>% 
        group_by(gene_id) %>%  
        summarise(rank=mean(value)) %>%
        # left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>%
        distinct() %>%
        arrange(-rank)
      rmats_res_gsea <- run_fgsea_analysis(ranks_rmats, is_ranked_already = F, organism = ORGANISM, gene_list_from = "rmats")
    }
    return(list(rmats_pos, background_genes))
}
 

combined_pos <- vector()
combined_neg <- vector()
null_background = F
for (tool in args$tool){
  
  if (tool == "vastools"){
    show("Processing vastools files..")
    genes <- process_vastools(args$vastools_positive_genes, args$vastools_negative_genes, args$vastools_ranks, args$vastools_ranks_groups, args$no_custom)
    combined_pos <- c(combined_pos, genes[[1]])
    if (is.null(args$vastools_negative_genes) || args$no_custom){
      null_background = T
    }
    else{
      combined_neg <- c(combined_neg, genes[[2]])
    }
    }
    
  else if (tool == "majiq"){
    show("Processing majiq files..")
    genes <- process_majiq(args$majiq_positive_genes, args$majiq_negative_genes, args$no_custom)
    combined_pos <- c(combined_pos, genes[[1]])
    if (is.null( args$majiq_negative_genes) || args$no_custom){
      null_background = T
    }
    else{
      combined_neg <- c(combined_neg, genes[[2]])
    }
    }
    
  else if (tool == "rmats"){
    show("Processing rmats files..")
    genes <- process_rmats(args$rmats_positive_genes, args$rmats_negative_genes, args$rmats_ranks, args$no_custom)
    combined_pos <- c(combined_pos, genes[[1]])
    if (is.null( args$majiq_negative_genes) || args$no_custom){
      null_background = T
    }
    else{
      combined_neg <- c(combined_neg, genes[[2]])
    }
  }
}

###########################
######TOOLS COMBINED ######
###########################
if (length(args$tool) > 1){
  show("Combining genes from multiple tools..")
  all_pos <- unique(combined_pos)
  show(paste0("Number of unique combined genes that will be used: ",length(all_pos)))
  if (null_background){
    all_neg <- NULL
  }
  else{
    all_neg <- unique(combined_neg)
  }
  all_res <- run_enrichment_gprofiler(all_pos, custom_genes_background = all_neg,
                                      retrieve_only_sign_results = T,
                                      exclude_iea = F,
                                      organism = ORGANISM,
                                      retrieve_short_link = F,
                                      measure_under = F, 
                                      domain_scope = "annotated")
  if(!is.null(all_res)){
    pdf(paste0(OUTBASENAME, "_gostplot_tools_combined.pdf"))
    p <- gostplot(all_res[[2]], interactive = F)
    plot(p)
    dev.off()
    
    pdf(paste0(OUTBASENAME, "_gprofiler_tools_combined.pdf"))
    p <- plot_enrichment_results(all_res[[1]], short_term_size=USE_SHORT_TERMS, size_of_short_term=SHORT_TERM_SIZE)
    plot(p)
    dev.off()
    write.table(all_res[[1]], file = paste0(OUTBASENAME, "_gprofiler_tools_combined.tsv"), row.names = F, quote=F, sep="\t")
  }
}
show("Done!")
