library(dplyr)
source("/Users/pbarbosa/git_repos/bioinfo_utils/r-scripts/rnaseq/standard_rna_seq.R")


#' @title GO enrichment of splicing events
#' @description Given a character vector of positive genes (with any relevant splicing event) 
#' performs gene list functional enrichment analysis using gprofiler or goseq packages.
#' 
#' If a set of background set of genes is provided, it will serve as a custom set of genes to 
#' fetch terms that will compose the statistical domain to perform the enrichment tests.
#' 
#' If a list of character vectors is provided in the \code{positive_genes} argument, a comparison
#' between gene lists will be performed by passing the \code{multi_query} argument to gprofiler.
#' It's useful to compare functional enrichment profiles between splicing events detected by
#' different tools.
#'  
#' The background list represents genes that meet similar read coverage criteria as the positive
#' gene list, so GO enrichment of highly expressed genes is avoided.
#'
#' @param positive_genes a character vector (or a list of character vectors) containing gene IDs that contain 
#' differentially spliced events.
#' @param background_genes a character vector containing gene IDs to serve as the background in the statistical
#' procedure (generated using vast-tools compare with --GO option set.). If not set, gprofiler will be used in 
#' the \code{package} argument because it can automatically fetch annotated lists of IDs for any given data
#' sources to enrich. As opposed to \code{positive_genes} this argument doesn't accept a list of several backrgound 
#' gene sets. Therefore, be careful on the choice of the background gene set when providing multiple positive gene lists.
#' @param package optional parameter that referes to the choice of the R package to run the tests. \code{gprofiler}
#' is used by default. Options: \code{gprofile, goseq}.
#' @param ... Arbitary number of additional arguments to be passed to the enrichment functions (check gprofiler and goseq 
#' enrichment functions to know which arguments can be passed). \code{run_enrichment_gprofiler} and \code{run_goseq_analysis}
#' 
#' @return a dataframe with the GO enriched terms
do_GO_functional_enrichment <- function(positive_genes, background_genes = NULL, package = c("gprofiler", "goseq"), ...) {
  
  package <- match.arg(package)
  if (is.null(background_genes) & package == "goseq"){
    stop("A background list of genes is required to run goseq package. Please set 'gprofiler' in the package argument
    or provide a list of genes in the 'backround_genes' argument.")
  } 
  else if (is.null(background_genes)){
    results <- run_enrichment_gprofiler(positive_genes, custom_genes_background = NULL, ...) 
  }
  else if (package == "gprofiler"){
    results <- run_enrichment_gprofiler(positive_genes, custom_genes_background = background_genes, ...) 
  }
  else if (package == "goseq"){
    results <- run_goseq_analysis(positive_genes, background_genes, "hg19")
  }
  return(results)
}

setwd("/Users/pbarbosa/MEOCloud/analysis/christian/rna_seq_analysis/human_CBE_mock/splicing/")
################################
########VAST-TOOLS##############
################################
#GO
vast_positive_genes <- read_tsv("vast-tools/2_positive_geneIDs_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()
vast_background_genes <- read_tsv("vast-tools/2_negative_geneIDs_to_GO.txt", col_names = "gene_id" ) %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

res <- do_GO_functional_enrichment(rmats_positive_genes, background_genes = NULL, package = "gprofiler",
                                          retrieve_only_sign_results = T,
                                          exclude_iea = F,
                                          retrieve_short_link = F,
                                          measure_under = F, 
                                          domain_scope = "annotated")


#FGSEA
vasttools_to_rank <- read_tsv("vast-tools/1_mock_CBE_TIDY_to_rank_genes.tsv") %>% separate(col = EVENT, into = c("symbol", "eventid"), sep="=")
vasttools_to_rank$mock <- rowMeans(subset(vasttools_to_rank, select = grepl("mock", colnames(vasttools_to_rank))), na.rm = T)
vasttools_to_rank$CBE <- rowMeans(subset(vasttools_to_rank, select = grepl("CBE", colnames(vasttools_to_rank))), na.rm = T)

ranks_vasttools <- vasttools_to_rank %>% mutate(dPSI = CBE - mock) %>%  group_by(symbol) %>% 
  drop_na(dPSI) %>%
  summarise(rank=max(dPSI)) %>%
  filter(symbol != "") %>%
  left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>%
  distinct() %>%
  arrange(-rank)

res <- run_fgsea_analysis(ranks_vasttools, is_ranked_already = F, gene_list_from = "vast-tools")


################################
########## MAJIQ ###############
################################
#GO
majiq_positive_genes <- read_tsv("majiq/mock_vs_CBE_majiq_positive_list_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()
majiq_background_genes <- read_tsv("majiq/mock_vs_CBE_majiq_negative_list_to_GO.txt", col_names = "gene_id" ) %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

res <- do_GO_functional_enrichment(majiq_positive_genes, background_genes = NULL, package = "gprofiler",
                                   retrieve_only_sign_results = T,
                                   exclude_iea = F, 
                                   retrieve_short_link = F,
                                   measure_under = F, 
                                   domain_scope = "annotated")
plot_enrichment_results(res)
#FGSEA
#Not robust at all. Ranks are made quite arbitarily as in majiq I select the maximum dPSI on a given junction
#from LSV. We can't really be sure of the directionality of the change.
majiq_to_rank <- read_tsv("majiq/mock_vs_CBE_majiq_gene_ranks_to_GSEA.txt")
ranks_majiq <- majiq_to_rank %>% mutate(rank=probability_max_dPSI * max_dPSI) %>%
  group_by(gene_id) %>% summarise(rank = mean(rank)) %>%
  left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>%
  distinct()

res <- run_fgsea_analysis(ranks_majiq, is_ranked_already = F, gene_list_from = "majiq")


################################
########## rMATS ###############
################################
#GO
rmats_positive_genes <- read_tsv("rmats/rmats_positive_list_to_GO.txt", col_names = "gene_id") %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()
rmats_background_genes <- read_tsv("rmats/rmats_negative_list_to_GO.txt", col_names = "gene_id" ) %>% 
  filter(str_detect(gene_id, "ENSG")) %>% pull()

res <- do_GO_functional_enrichment(rmats_positive_genes, background_genes = NULL, package = "gprofiler",
                                   retrieve_only_sign_results = T,
                                   exclude_iea = F,
                                   retrieve_short_link = F,
                                   measure_under = F, 
                                   domain_scope = "annotated")
plot_enrichment_results(res)

#FGSEA
rmats_to_rank <- read_tsv("rmats/rmats_gene_ranks_to_GSEA.txt", col_names = c("gene_id", "padj", "dPSI"))
ranks_rmats <- rmats_to_rank %>% mutate(padj = ifelse(padj == 1, 0.9999999, ifelse(padj == 0, 1e-20, padj)),
                      dPSI = -dPSI, value = -log10(padj) * dPSI) %>%  group_by(gene_id) %>% 
                      summarise(rank=mean(value)) %>%
                      left_join(dplyr::select(ensembl_genes, c(gene_id, symbol))) %>%
                      distinct() %>%
                      arrange(-rank)

res <- run_fgsea_analysis(ranks_rmats, is_ranked_already = F, npermutations = 2000, gene_list_from = "rmats")



###############################
#######TOOLS COMBINED #########
###############################
all_positive_genes <- Reduce(dplyr::union, list(vast_positive_genes, majiq_positive_genes, rmats_positive_genes)) 
all_background_genes <- Reduce(dplyr::union, list(vast_background_genes, majiq_background_genes, rmats_background_genes))

res <- do_GO_functional_enrichment(all_positive_genes, background_genes = NULL, package = "gprofiler",
                                   retrieve_only_sign_results = T,
                                   exclude_iea = F,
                                   retrieve_short_link = F,
                                   measure_under = F, 
                                   domain_scope = "annotated")

plot_enrichment_results(res, short_term_size = F, top_n = 20, size_of_short_term = 500, add_info_to_labels = T)
dim(res)
######################
####COMPARE TOOLS ####
######################
res <- do_GO_functional_enrichment(list(rmats_positive_genes, majiq_positive_genes, vast_positive_genes), background_genes = NULL, package="gprofiler",
                                   multi_query = T,
                                   retrieve_only_sign_results = T, 
                                   exclude_iea = F, 
                                   retrieve_short_link = F,
                                   measure_under = F, 
                                   domain_scope = "annotated"
                                  )


to_heatmap <- res %>% dplyr::select(c(term_name, significant)) %>%
  column_to_rownames("term_name") %>%
  separate(col = significant, sep="\\,", into = c("rMATS", "MAJIQ", "VAST-TOOLS"), convert = T) 
  #mutate_if(is.character, funs(. == "FALSE"))
heatmap(test, scale = "none", Rowv = NA, Colv = NA, col = cm.colors(2), main = "HeatMap Example") 
