library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(igraph)
library(ggraph)
library(tidyverse)
# library(ViSEAGO)
library(tibble)
library(dplyr)
# library(dbplyr)
library(reshape2)
library(DOSE)

# wrapper functions -------------------

#' Wrapper function which perform GSEA and over-representation analyses using
#' ClusterProfiler for the dataset submitted.
#' 
#' Wrapper function which produces a GSEA analysis and two over-representation analyses
#' The two over-representation analyses are done on down-regulated and up-regulated genes
#' respectively.
#' 
#' The function then calls a plot-wrapper function which produces all the plots for
#' the GSEA and over-representation analyses.
#'
#' @param DE_file The differential expression results
#' @param save_path The path used to save the results
#' @return None
#'
#' @examples
#' DE_file <- read.csv("DE_results/PBMC_2_vs_IgD/DE_results_with_gtf.csv"),check.names = FALSE)
#' save_path=paste0("cluster_profiler_results/PBMC_2_vs_IgD")
#' wrapper_for_GSEA_and_overrepresentation(DE_file,save_path)
#'
#' @export
wrapper_for_GSEA_and_overrepresentation <- function(DE_file,filter_choice,save_path,selected_org_gsea,selected_org_GO,go_dta=NULL,exp_name=''){
  
  
  gse <- GSEA_wrapper(DE_file, filter_on =filter_choice, selected_org_gsea,FALSE, key_type = "SYMBOL")
  go_enrich_lst <- GSEA_wrapper(DE_file, filter_on =filter_choice,selected_org_gsea,TRUE,ontology = "BP", key_type = "SYMBOL")
  go_enrich_up <- go_enrich_lst[[1]]
  go_enrich_down <- go_enrich_lst[[2]]
  gene_list <- DE_file$log2FoldChange
  names(gene_list) <- DE_file$gene_id
  gene_list <- na.omit(gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  if (is.null(go_dta)==T){
    d <- godata(selected_org_GO, ont="BP")
  }else{
    d<-go_dta
  }

  gse <- pairwise_termsim(gse, method="Wang", semData = d)
  results_over_up <- safe_over_rep(go_enrich_up,d,save_path,'up')
  go_enrich_up <- results_over_up[[1]]
  enrich_up <- results_over_up[[2]]
  
  results_over_down <- safe_over_rep(go_enrich_down,d,save_path,'down')
  go_enrich_down <- results_over_down[[1]]
  enrich_down <- results_over_down[[2]]
  
  custom_plots(gse,organism_str = selected_org_GO,organism_obj = selected_org_gsea,save_path=save_path,ontology = 'BP',exp_name=exp_name)
  
  if (enrich_up==T & enrich_down==T){
    plot_wrapper_GSEA_over_rep(gse,go_enrich_up,go_enrich_down,gene_list,save_path)
    save(list = c("gse","go_enrich_up","go_enrich_down","gene_list"), file = paste0(save_path,"/cluster_profiler_objects.RData"))
  }else if (enrich_up==T & enrich_down==F){
    plot_wrapper_GSEA_over_rep(gse,go_enrich_up,go_enrich_down=NULL,gene_list,save_path)
    save(list = c("gse","go_enrich_up","gene_list"), file = paste0(save_path,"/cluster_profiler_objects.RData"))
  }else if (enrich_up==F & enrich_down==T){
    plot_wrapper_GSEA_over_rep(gse,go_enrich_up=NULL,go_enrich_down,gene_list,save_path)
    save(list = c("gse","go_enrich_down","gene_list"), file = paste0(save_path,"/cluster_profiler_objects.RData"))
  }else{
    plot_wrapper_GSEA_over_rep(gse,go_enrich_up=NULL,go_enrich_down=NULL,gene_list,save_path)
    save(list = c("gse","gene_list"), file = paste0(save_path,"/cluster_profiler_objects.RData"))
  }
  

  # save.image(paste0(save_path,"/cluster_profiler_objects.RData"))#Saves all objects in memory in a binary file
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Create a standard or overexpression GSEA analysis object
#' 
#' The function takes in a dataframe containing differential expression results 
#' Using the organism, ontology, and key_type provided, the function
#' will calculate and evaluate the significance of GO terms found. A pvalue cutoff can be given 
#' to customize significant GO terms.
#'
#' The function will perform a standard or overexpression analysis based on the users
#' choice through the overexpression parameter
#'
#' @param DE_file Differential expression results
#' @param organism The annotation object for a specific organism
#' @param overexpression A boolean indicating if the analysis will be standard or overexpression
#' @param ontology The type of GO term to be found, default is ALL
#' @param key_type The source of the annotation, typing should match gene_id typing of DE_file
#' @param pval_cutoff The significance cutoff for pvalue
#' @param qval_cutoff The significance cutoff for qval
#' 
#' @return gseaResults object
#'
#' @examples
#' DE_res <- read.csv("DE_results/PBMC_2_vs_IgD/DE_results_with_gtf.csv"),check.names = FALSE)
#' gse <- standard_GSEA_wrapper(DE_res, org.Dm.eg.db, FALSE)
#'
#' @export
GSEA_wrapper <- function (DE_file, filter_on,organism, overexpression, ontology="BP",key_type="SYMBOL", FDR_cutoff=0.05){
  original_gene_list <- DE_file$log2FoldChange
  names(original_gene_list) <- DE_file$gene_id
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(original_gene_list, decreasing = TRUE)
  
  my_universe <- keys(organism,keytype = key_type)
  
  
  if (overexpression==TRUE){
    sig_genes_df <- DE_file[DE_file[filter_on] < FDR_cutoff,]
    # sig_genes_df = subset(DE_file, pvalue < pval_cutoff)
    pval_cutoff <- max(sig_genes_df$pvalue)
    
    genes <- sig_genes_df$log2FoldChange
    names(genes) <- sig_genes_df$gene_id
    
    genes <- na.omit(genes)
    
    genes_up <- names(genes)[genes > 0]
    genes_down <- names(genes)[genes < 0]
    

    gse_up <- enrichGO(gene = genes_up,
                       universe = my_universe,
                       OrgDb = organism, 
                       keyType = key_type,
                       readable = F,
                       ont = ontology,
                       pvalueCutoff = pval_cutoff,
                       qvalueCutoff = FDR_cutoff)
    
    gse_down <- enrichGO(gene = genes_down,
                         universe = my_universe,
                         OrgDb = organism, 
                         keyType = key_type,
                         readable = F,
                         ont = ontology,
                         pvalueCutoff = pval_cutoff,
                         qvalueCutoff = FDR_cutoff)
    gse <- list(gse_up,gse_down)
  }else{
    gse <- gseGO(geneList=gene_list, 
                 ont =ontology, #Can be BP (biological process), MF (molecular function), CC(Cellular Component) or all
                 keyType = key_type, #Type of annotation used, many possibilities
                 minGSSize = 3, #minimum number of genes in the set
                 maxGSSize = 800, #maximum number of genes in the set
                 pvalueCutoff = 0.05, #pvalue cutoff
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none") #several possibilities
    # if (dim(gse)[1]==0){
    #   print("No enriched terms were found, possibly due to a strict threshold")
    #   return()
    # }else{
    #   return (gse)
    # }
    
  }
}
# Calculation functions -------------------

#' A wrapper function to check if there are any over-representation results
#' 
#' The GSEA/ClusterProfiler pipeline produces go_enrich, which contains the 
#' over-representation results from cluster profiler. This function checks if the 
#' the go_enrich object can be used via the calculation of pairwise_termsmim.
#' If the function goes through, the necessary folders are created and a boolean
#' for True is returned along with the pairwise results.
#' If the function does not succeed, a False boolean is returned along with an empty
#' list.
#'
#' @param go_enrich The overrepresentation results as produced by cluster profiler
#' @param d The results from the \code{godata} function
#' @param save_path The path where the folders will be created
#' @param type Either 'up' or 'down' indicating if the over-representation
#' was done with the up-regulated genes or the down-regulated genes
#' 
#' @return results A list containing the go_enrich results along with a boolean
#' stating if the results are empty (FALSE) or not (TRUE)
#'
#' @examples
#' DE_res <- read.csv("DE_results/PBMC_2_vs_IgD/DE_results_with_gtf.csv"),check.names = FALSE)
#' gse <- standard_GSEA_wrapper(DE_res, org.Dm.eg.db, FALSE)
#'
#' @export
safe_over_rep <- function(go_enrich, d,save_path,type){
  normal_over_rep <- function(go_enrich, d,save_path,type){
    go_enrich_res <- pairwise_termsim(go_enrich, method="Wang", semData = d)
    enrich_bool<-T
    
    if (type=='up'){
      dir.create(paste0(save_path,"/Over_representation_up"))
    }else if (type=='down'){
      dir.create(paste0(save_path,"/Over_representation_down"))
    }
    return(list(go_enrich_res,enrich_bool))
  }
  bad_over_rep <- function(err){
    go_enrich_res <- NULL
    enrich_bool <- F
    return(list(go_enrich_res,enrich_bool))
  }
  
  results <- tryCatch(normal_over_rep(go_enrich, d,save_path,type), error=bad_over_rep)
  return(results)
}






# Plotting functions -------------------

#' A series of png calls to produce and save all the plots required
#' 
#' This wrapper function takes in one GSEA and two over-representation results
#' and plots all the necessary plots for each analysis. The pots are stored in 
#' folders which represent their analysis.
#' 
#' @param gse The GSEA object as produced by ClusterProfiler's \code{gseGO} function
#' @param go_enrich_up The over-representation object for the up-regulated significant
#' genes as produced by ClusterProfiler's \code{GOenrich} function
#' @param go_enrich_down The over-representation object for the down-regulated significant
#' genes as produced by ClusterProfiler's \code{GOenrich} function
#' @param gene_list a gene list to be used to create the fold change legend in cnetplots
#' 
#' @return None
#' 
#' @example 
#' DE_file <- read.csv("DE_results/PBMC_2_vs_IgD/DE_results_with_gtf.csv",check.names = FALSE)
#' save_path=paste0("cluster_profiler_results/PBMC_2_vs_IgD")
#' 
#' gse <- GSEA_wrapper(DE_file, org.Hs.eg.db,FALSE, key_type = "SYMBOL")
#'   
#' go_enrich_lst <- GSEA_wrapper(DE_file, org.Hs.eg.db,TRUE,ontology = "BP", key_type = "SYMBOL")
#' go_enrich_up <- go_enrich_lst[[1]]
#' go_enrich_down <- go_enrich_lst[[2]]
#' 
#' gene_list <- DE_file$log2FoldChange
#' names(gene_list) <- DE_file$gene_id
#' gene_list <- na.omit(gene_list)
#' gene_list = sort(gene_list, decreasing = TRUE)
#' 
#' plot_wrapper_GSEA_over_rep(gse,go_enrich_up,go_enrich_down,gene_list)
#' 
#' @export
plot_wrapper_GSEA_over_rep <- function(gse,go_enrich_up,go_enrich_down,gene_list,save_path){

  
  
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/dotplot.png"),width = 30, height = 40,
  #                     plot=dotplot(gse, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign))
  # 
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/emapplot_condensed.png"),
  #                     plot=emapplot(gse, showCategory = 50,color="enrichmentScore"))
  # 
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/emapplot_enhanced.png"),width=50,height=50,
  #                       plot=emapplot(gse, showCategory = 50,color="enrichmentScore"))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/cnetplot.png"),
                      plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = 1))
  
  
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/gseaplot.png"),
  #                     plot=gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1))
  # 
  # 
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/ridgeplot.png"),width = 30, height = 40,
  #                     plot=ridgeplot(gse, showCategory = 10) + labs(x = "enrichment distribution"))
  # 
  # my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/pmcplot.png"),
  #                     plot=pmcplot(gse$Description[1:3], 2010:2019, proportion=FALSE))
  
  if (is.null(go_enrich_up)==F){
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/barplot.png"),
                        plot=barplot(go_enrich_up, drop = TRUE, showCategory = 20, title = "GO Biological Pathways",font.size = 8))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/dotplot.png"),width = 30, height = 40,
                        plot=dotplot(go_enrich_up, showCategory = 20, title = "Enriched Pathways"))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/emapplot_condensed.png"),
                        plot=emapplot(go_enrich_up, showCategory = 50))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/emapplot_enhanced.png"),width=50,height=50,
                        plot=emapplot(go_enrich_up, showCategory = 50))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/cnetplot_condensed.png"),
                        plot=cnetplot(go_enrich_up, categorySize="pvalue", foldChange=gene_list, showCategory = 5))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/cnetplot_enhanced.png"),width=50,height=50,
                        plot=cnetplot(go_enrich_up, categorySize="pvalue", foldChange=gene_list, showCategory = 5))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/goplot_condensed.png"),
                        plot=goplot(go_enrich_up, showCategory = 10))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_up/goplot_enhanced.png"),width = 75, height=75,
                        plot=goplot(go_enrich_up, showCategory = 10))

  }
  
  if (is.null(go_enrich_down)==F){
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/barplot.png"),
                        plot=barplot(go_enrich_down, drop = TRUE, showCategory = 20, title = "GO Biological Pathways",font.size = 8))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/dotplot.png"),width = 30, height = 40,
                        plot=dotplot(go_enrich_down, showCategory = 20, title = "Enriched Pathways"))
    
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/emapplot_condensed.png"),
                        plot=emapplot(go_enrich_down, showCategory = 50))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/emapplot_enhanced.png"),width=50,height=50,
                        plot=emapplot(go_enrich_down, showCategory = 50))
    
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/cnetplot_condensed.png"),
                        plot=cnetplot(go_enrich_down, categorySize="pvalue", foldChange=gene_list, showCategory = 5))
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/cnetplot_enhanced.png"),width=50,height=50,
                        plot=cnetplot(go_enrich_down, categorySize="pvalue", foldChange=gene_list, showCategory = 5))
    

    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/goplot_condensed.png"),
                        plot=goplot(go_enrich_down, showCategory = 10))
    
    
    my_gg_save_function(save_path_name = paste0(save_path,"/Over_representation_down/goplot_enhanced.png"),width = 75, height=75,
                        plot=goplot(go_enrich_down, showCategory = 10))
    

  }

}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A modified function from enrichplot's goplot function
#' 
#' Function is modified to allow the use of this function for gsea analyses (previously
#' only over/under representations were accepted) It also enables the ability to submit 
#' GO IDs as an input instead of a number which reflects the top n pathways.
#' 
#' @param x The GSEA or GOEnrich object
#' @param showCategory a number to plot the top n pathways
#' @param color adapted to take in either "p.adjust" or "enrichmentScore"
#' @param layout layout used during the plot as per ggplot2
#' @param geom graphical input for geom
#' @param selected_IDs a list of GO IDs to be plotted, if not null, it will override 'showCategory'
#' 
#' @return None
#' 
#' @example 
#' goplot_modified(
#' gse,
#' selected_IDs=GO_list,
#' color = "enrichmentScore"
#' )
#' 
#' 
#' @export
goplot_modified <- function(x, showCategory = 10, color = "p.adjust",layout = "sugiyama", geom = "text", selected_IDs=NULL) {
  GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
  getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
  if (is.null(selected_IDs)==T){
    # Replaced function call with function content
    n <- showCategory
    if (nrow(x) < n) {
      n <- nrow(x)
    }
    
    y <- as.data.frame(x)
    y <- y[1:n,]
    id <- y$ID[1:n]
  }else{
    id <- selected_IDs
  }
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  
  if (!exists(".GOSemSimEnv")) GOSemSim_initial()
  .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
  gotbl <- get("gotbl", envir=.GOSemSimEnv)
  
  GOANCESTOR <- getAncestors(x@ontology)
  anc <- AnnotationDbi::mget(id, GOANCESTOR)
  ca <- anc[[1]]
  for (i in 2:length(anc)) {
    ca <- intersect(ca, anc[[i]])
  }
  
  uanc <- unique(unlist(anc))
  uanc <- uanc[!uanc %in% ca]
  dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]
  
  
  edge <- dag[, c(5, 1, 4)]
  node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
  node$color <- x[node$go_id, color]
  node2 <- node

  node$size <- sapply(geneSets[node$go_id], length)
  g <- graph.data.frame(edge, directed=TRUE, vertices=node)
  
  E(g)$relationship <- edge[,3]
  
  p <- ggraph(g, layout=layout) +
    ## geom_edge_link(aes_(color = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm')) +
    geom_edge_link(aes_(linetype = ~relationship),
                   arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm'),
                   colour="darkgrey") +
    ## geom_node_point(size = 5, aes_(fill=~color), shape=21) +
    geom_node_point(size = 5, aes_(color=~color)) +
    theme_void() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=TRUE))
  ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE))
  
  if (geom == "label") {
    p <- p + geom_node_label(aes_(label=~Term, fill=~color), repel=TRUE) +
      scale_fill_continuous(low="red", high="blue", name = color,
                            guide=guide_colorbar(reverse=TRUE), na.value="white")
    ## scale_fill_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
  } else {
    p <- p + geom_node_text(aes_(label=~Term), repel=TRUE)
  }
  return(p)
}


# Custom plotting functions -------------------

#' Function to save non-gg based plots using ggsave function
#' 
#' This function was developed to utilize the superior ggsave function
#' as opposed to the sub-optimal png function. It uses a try catch as non gg plots
#' cause an grid.draw error which has NO impact, but halts the script.
#' 
#' The function can optionally take in a custom width and height.
#' 
#' @param save_path_name The path (also containing the file name) in which
#' the plot is to be saved
#' @param plot the plot
#' @param height optional height in cm
#' @param width optional width in cm
#' 
#' @return None
#' 
#' @example 
#' 
#' 
#' 
#' @export
my_gg_save_function <- function(save_path_name,plot,height=NULL,width=NULL){
  tryCatch({
    if (is.null(height)==F){
      ggsave(save_path_name,dpi=300,height=height,width=width,unit='cm',
             plot,bg = "white")
    }else{
      ggsave(save_path_name,dpi=300,
             plot,bg = "white")
    }
    
  },
  error=function(cond){print(cond)}
  )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Wrapper function to save cnet and emaps using custom gsea calls
#' 
#' Simple function which calls my_gg_save to create custom emap and cnet plots
#' based on the descriptions within GSEA results.
#' 
#' @param gse The GSEA object as resulting from ClusterPofiler
#' @param desc_list a vector containing the descriptions to be plotted
#' @param save_path The path in which each plot will be saved
#' 
#' @return None
#' 
#' @example 
#' 
#' 
#' 
#' @export
plot_cnet_emap_go_using_gsea_desc <- function(gse,desc_list,save_path){
  my_gg_save_function(paste0(save_path,'/custom_emap.png'),plot=emapplot(gse,showCategory =my_desc_list))
  
  GO_list <- gse@result$ID[gse@result$Description %in% my_desc_list]
  attr(gse,"ontology") <- "BP"
  my_gg_save_function(paste0(save_path,'/','custom_goplot_enrich.png'),
                      plot=goplot_modified(
                        gse,
                        selected_IDs=GO_list,
                        # color = "p.adjust",
                        color = "enrichmentScore",
                        layout = "sugiyama",
                        geom = "text"
                      ),height=50,width=30)
  
  my_gg_save_function(paste0(save_path,'/','custom_goplot_significance.png'),
                      plot=goplot_modified(
                        gse,
                        selected_IDs=GO_list,
                        color = "p.adjust",
                        # color = "enrichmentScore",
                        layout = "sugiyama",
                        geom = "text"
                      ),height=50,width=30)
  
  for (desc in my_desc_list){
    my_gg_save_function(paste0(save_path,'/',desc,'_cnetplot.png'),plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = desc, node_label='gene'))
  }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A modified version of cnetplot from the enrichplot package
#' 
#' 
#' @param x the geneset
#' @param showcategory a numerical or vector of descriptions to plot
#' @param foldchange a named vector with genes and associated logfoldchange
#' @param colorEdge Logical, whether coloring edge by enriched terms, the default value is FALSE. 
#' @param circular Logical, whether using circular layout, the default value is FALSE.
#' @param node_label Select which labels to be displayed.
#' one of 'category', 'gene', 'all'(the default) and 'none'.
#' @param cex_category Number indicating the amount by which plotting category
#' nodes should be scaled relative to the default, the default value is 1.
#' @param cex_gene Number indicating the amount by which plotting gene nodes
#' should be scaled relative to the default, the default value is 1.
#' @param cex_label_category Scale of category node label size, the 
#' default value is 1.
#' @param cex_label_gene Scale of gene node label size, the default value is 1.
#' 
#' @return 
#' 
#' @example 
#' 
#' 
#' 
#' @export
cnetplot_mod <- function(x,
                         showCategory = 5,
                         foldChange   = NULL,
                         layout = "kk",
                         colorEdge = FALSE,
                         circular = FALSE,
                         node_label = "all",
                         cex_category = 1,
                         cex_gene = 1,
                         node_label_size = NULL,
                         cex_label_category = 1,
                         cex_label_gene = 1,
                         ...) {
  
  if (!is.null(node_label_size))
    message("node_label_size parameter has been changed to 'cex_label_category' and 'cex_label_gene'")
  
  label_category <- 5
  label_gene <- 5
  node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  } else {
    geom_edge <- geom_edge_link
  }
  
  # geneSets <- extract_geneSets(x, showCategory)
  geneSets <- x
  
  g <- list2graph_mod(geneSets)
  
  # foldChange <- fc_readable(x, foldChange)
  
  
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
  } else {
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
  }
  
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(n+1):length(V(g))] <- fc
    show_legend <- c(TRUE, FALSE)
    names(show_legend) <- c("color", "size")
    p <- ggraph(g, layout=layout, circular = circular)
    p <- p + edge_layer +
      # geom_node_point(aes_(color=~as.numeric(as.character(color)),
      geom_node_point(aes_(color=~I("#E5C494"), size=~size),
                      data = p$data[1:n, ]) +
      scale_size(range=c(3, 8) * cex_category) +
      ggnewscale::new_scale("size") +
      ggnewscale::new_scale_color() +
      geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                      data = p$data[-(1:n), ], show.legend = show_legend) +
      scale_size(range=c(3, 3) * cex_gene) +
      scale_colour_gradient2(name = "fold change", low = "blue",
                             mid = "white", high = "red")
  } else {
    V(g)$color <- "#B3B3B3"
    V(g)$color[1:n] <- "#E5C494"
    p <- ggraph(g, layout=layout, circular=circular)
    p <- p + edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size), data = p$data[1:n, ]) +
      scale_size(range=c(3, 8) * cex_category) +
      ggnewscale::new_scale("size") +
      geom_node_point(aes_(color=~I(color), size=~size),
                      data = p$data[-(1:n), ], show.legend = FALSE) +
      scale_size(range=c(3, 3) * cex_gene)
  }
  
  p <- p + theme_void()
  
  # Unlikely to have to touch this
  if (node_label == "category") {
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                              size = label_category * cex_label_category, bg.color = "white")
    } else {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,],
                              size = label_category * cex_label_category)
    }
  } else if (node_label == "gene") {
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                              repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white")
    } else {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                              repel=TRUE, size = label_gene * cex_label_gene)
    }
  } else if (node_label == "all") {
    if (utils::packageVersion("ggrepel") >= "0.9.0") {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                              repel=TRUE, size = label_gene * cex_label_gene, bg.color = "white") + 
        geom_node_text(aes_(label=~name), repel=TRUE,
                       size = label_category * cex_label_category, bg.color = "white", data = p$data[1:n,]) 
    } else {
      p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                              repel=TRUE, size = label_gene * cex_label_gene) + 
        geom_node_text(aes_(label=~name), data = p$data[1:n,],
                       repel=TRUE, size = label_category * cex_label_category)
    }
    
  }
  
  return(p)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A modified version of emapplot from the enrichplot package to allow it's
#' use for the processed version of the visEAGO results
#' 
#' The function previously contained several elements which check if the object
#' was a clusterProfiler object, it also created the GeneSets and similarity matrices
#' form this object. These elements were removed and replced with custom inputs of genesets
#' and similarity matrices.
#' 
#' @param y_custom A dataframe of GO pathway results
#' @param showCategory A numerical value or list of descriptions to plot
#' @param simTerm_table A similarity matrix for the GO pathways to be plotted
#' @param geneSets A geneset (all pathways and genes within those pathways)
#' @param color a string, the column name of y for nodes colors
#' @param layout the type of layout to use
#' @param min_edge minimum percentage of overlap genes to display the edge,
#' should between 0 and 1, default value is 0.2
#' @param cex_label_category size of node label
#' @param cex_label_category scale of category node label size
#' @param cex_category number indicating the amount by which plotting category
#' @param shadowtext if the text should be shadowed
#' 
#' @return 
#' 
#' @example 
#' 
#' 
#' 
#' @export
emapplot_visEAGO <- function(y_custom, showCategory = 30,simTerm_table,geneSets, color="pvalue",
                             layout = "nicely", min_edge=0.2,
                             cex_label_category  = 1, cex_category = 1,
                             cex_line = 1, shadowtext = TRUE) {
  
  label_category <- 5
  n <- showCategory
  
  y <- as.data.frame(y_custom)
  g <- get_igraph_mod(simTerm_table=simTerm_table,geneSets=geneSets, y=y, n=n, color=color, cex_line=cex_line,
                      min_edge=min_edge)
  
  # if(n == 1) {
  #   return(ggraph(g,"tree") + geom_node_point(color="red", size=5) +
  #            geom_node_text(aes_(label=~name)))
  # }
  
  p <- ggraph(g, layout=layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
  }
  
  p <- ggraph(g, layout=layout)
  if (length(E(g)$width) > 0) {
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
  }
  p <- p + geom_node_point(aes_(color=~color, size=~size))
  
  if (utils::packageVersion("ggrepel") >= "0.9.0") {
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                            size = label_category * cex_label_category, bg.color = "white")
  } else {
    p <- p + geom_node_text(aes_(label=~name), repel=TRUE,
                            size = label_category * cex_label_category)
  }
  # geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
  if (color=='enrichmentScore'){
    p + theme_void() +
      scale_fill_gradientn(colours=rev(c("red","white","blue")),
                           # values=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                           breaks=c(-1,0,1),
                           # breaks=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                           limits=c(-1, 1), oob=squish,
                           guide=guide_colorbar(reverse=F), name = color,aesthetics = "colour")+ 
      scale_size(range=c(3, 8) * cex_category)
  }else{
    p + theme_void() +
      scale_color_continuous(low="red", high="blue", name = color,
                             # breaks=c(min(V(g)$color),max(V(g)$color)),
                             # label=c(min(V(g)$color),max(V(g)$color)),
                             guide=guide_colorbar(label=T,reverse=F)) +
      scale_size(range=c(3, 8) * cex_category)
  }
  
  

  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' A modified version of cnetplot from the enrichplot package
#' 
#' @param ontology the ontology
#' @param my_y GO pathways results
#' @param my_gene_sets The geneset
#' @param showCategory a numerical or vector of descriptions to plot
#' @param color a string, the column name of y for nodes colors
#' @param layout The layout to use
#' @param geom either geom text or geom label
#' @param selected_IDs The specific IDs to plot
#' 
#' @return 
#' 
#' @example 
#' 
#' 
#' 
#' @export
goplot_modified_visEAGO <- function(ontology, my_y, my_gene_sets, showCategory = 10, color = "p.adjust",layout = "sugiyama", geom = "text", selected_IDs=NULL) {
  GOSemSim_initial <- getFromNamespace(".initial", "GOSemSim")
  getAncestors <- getFromNamespace("getAncestors", "GOSemSim")
  
  id <- selected_IDs
  
  geneSets <- my_gene_sets ## use core gene for gsea result
  
  if (!exists(".GOSemSimEnv")) GOSemSim_initial()
  .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
  gotbl <- get("gotbl", envir=.GOSemSimEnv)
  
  GOANCESTOR <- getAncestors(ontology)
  anc <- AnnotationDbi::mget(id, GOANCESTOR)
  ca <- anc[[1]]
  for (i in 2:length(anc)) {
    ca <- intersect(ca, anc[[i]])
  }
  
  uanc <- unique(unlist(anc))
  uanc <- uanc[!uanc %in% ca]
  dag <- gotbl[gotbl$go_id %in% unique(c(id, uanc)),]
  
  
  edge <- dag[, c(5, 1, 4)]
  node <- unique(gotbl[gotbl$go_id %in% unique(c(edge[,1], edge[,2])), 1:3])
  
  node$color <- my_y[node$go_id, color]
  
  node$size <- sapply(geneSets[node$go_id], length)
  
  
  
  g <- graph.data.frame(edge, directed=TRUE, vertices=node)
  
  E(g)$relationship <- edge[,3]
  
  p <- ggraph(g, layout=layout) +
    ## geom_edge_link(aes_(color = ~relationship), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm')) +
    geom_edge_link(aes_(linetype = ~relationship),
                   arrow = arrow(length = unit(2, 'mm')), end_cap = circle(2, 'mm'),
                   colour="darkgrey") +
    ## geom_node_point(size = 5, aes_(fill=~color), shape=21) +
    geom_node_point(size = 5, aes_(color=~color)) +
    theme_void() +
    scale_color_continuous(low="red", high="blue", name = color,
                           # breaks=c(min(V(g)$color),max(V(g)$color)),
                           # label=c(min(V(g)$color),max(V(g)$color)),
                           guide=guide_colorbar(reverse=TRUE))
  
  
  ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE))
  
  if (geom == "label") {
    p <- p + geom_node_label(aes_(label=~Term, fill=~color), repel=TRUE) +
      scale_fill_continuous(low="red", high="blue", name = color,
                            guide=guide_colorbar(reverse=TRUE), na.value="white")
    ## scale_fill_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE), na.value='white')
  } else {
    p <- p + geom_node_text(aes_(label=~Term), repel=TRUE)
  }
  return(p)
}

library(scales)
dotplot_YL_mod_enrich_legend <- function(object, x = "geneRatio", color = "p.adjust",
                                 showCategory=10, size=NULL, split = NULL,
                                 font.size=12, title = "", orderBy="x",
                                 label_format = 50, decreasing=TRUE) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue","enrichmentScore"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split=split)
  
  #Replcae activated and suppressed by better terminology
  df$.sign<-gsub(x=df$.sign,pattern = 'activated',replacement = 'positively enriched')
  df$.sign<-gsub(x=df$.sign,pattern = 'suppressed',replacement = 'negatively enriched')
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  # #Set-up necessary functions
  # default_labeller <- function(n) {
  #   function(str){
  #     str <- gsub("_", " ", str)
  #     ep_str_wrap(str, n)
  #   }
  # }
  # 
  # ep_str_wrap <- function(string, width) {
  #   x <- gregexpr(' ', string)
  #   vapply(seq_along(x),
  #          FUN = function(i) {
  #            y <- x[[i]]
  #            n <- nchar(string[i])
  #            len <- (c(y,n) - c(0, y)) ## length + 1
  #            idx <- len > width
  #            j <- which(!idx)
  #            if (length(j) && max(j) == length(len)) {
  #              j <- j[-length(j)]
  #            }
  #            if (length(j)) {
  #              idx[j] <- len[j] + len[j+1] > width
  #            }
  #            idx <- idx[-length(idx)] ## length - 1
  #            start <- c(1, y[idx] + 1)
  #            end <- c(y[idx] - 1, n)
  #            words <- substring(string[i], start, end)
  #            paste0(words, collapse="\n")
  #          },
  #          FUN.VALUE = character(1)
  #   )
  # }
  # 
  
  
  # label_func <- default_labeller(label_format)
  # if(is.function(label_format)) {
  #   label_func <- label_format
  # }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  # df$Description<-gsub('(.{1,60})(\\s|$)', '\\1\n', df$Description)
  # df$Description<-paste(strwrap(df$Description,30), collapse="\n")
  new_col<-c()
  for(item in df$Description){
    new_item<-paste(strwrap(item,60), collapse="\n")
    new_col<-c(new_col,new_item)
  }
  df$Description<-new_col
  
  df$Description <- factor(df$Description,
                           levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_fill_gradientn(colours=rev(c("red","white","blue")),
                         # values=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                         breaks=c(-1,0,1),
                         # breaks=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                         limits=c(-1, 1), oob=squish,
                         guide=guide_colorbar(reverse=F), name = color,aesthetics = "colour")+ 
    # scale_color_continuous(low="blue", high="gray", name = color,
    #                        limits=c(-1, 0),
    #                        guide='legend') +
    scale_y_discrete(labels = function(x) str_wrap(x, width = label_format))+
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
    scale_size(range=c(3, 8))
  
}




#' A wrapper function to plot overrepresentation results of WGCNA modules
#' 
#' @param gse_over The overrepresentation object
#' @param save_path The path to which the plots are saved
#' 
#' @example 
#' 
#' 
#' 
#' @export
plot_overrep_WGCNA <- function(gse_over,save_path){
  my_gg_save_function(save_path_name = paste0(save_path,"/barplot.png"),
                      plot=barplot(gse_over, drop = TRUE, showCategory = 20, title = "GO Biological Pathways",font.size = 8))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/dotplot.png"),width = 30, height = 40,
                      plot=dotplot(gse_over, showCategory = 20, title = "Enriched Pathways"))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/emapplot_condensed.png"),
                      plot=emapplot(gse_over, showCategory = 50))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/emapplot_enhanced.png"),width=50,height=50,
                      plot=emapplot(gse_over, showCategory = 50))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/cnetplot_condensed.png"),
                      plot=cnetplot(gse_over, categorySize="pvalue", showCategory = 5))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/cnetplot_enhanced.png"),width=50,height=50,
                      plot=cnetplot(gse_over, categorySize="pvalue", showCategory = 5))
  
  
  my_gg_save_function(save_path_name = paste0(save_path,"/goplot_condensed.png"),
                      plot=goplot(gse_over, showCategory = 10))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/goplot_enhanced.png"),width = 75, height=75,
                      plot=goplot(gse_over, showCategory = 10))
}


#' A wrapper function to create a set of custom GSEA/overrepresentation plots
#' 
#' The function serves as a method to save enrichment maps and dotplots that use
#' enrichment score as a legend instead of the standard pvalue.
#' 
#' The plots can only be generated by using GSEA data as it relies on the presence
#' of enrichment score
#'
#' @param gse The GSEA object as produced by ClusterProfiler's \code{gseGO} function
#' @param organism_str A string of the annotation dbi organism to be used
#' @param organism_obj The annotation dbi object ot be used.
#' @param save_path The path to which the plots will be saved
#' @param ontology The ontology to be used
#' @param exp_name The title that will be given to the dotplots
#' 
#' @return None
#'
#' @examples
#'
#' @export
custom_plots <-function(gse,organism_str,organism_obj,save_path,ontology='BP',exp_name=''){
  
  Bioconductor<-ViSEAGO::Bioconductor2GO()
  myGENE2GO<-ViSEAGO::annotate(
    organism_str,
    Bioconductor)
  
  #Set up for dotplots
  gse_neg<-gse@result$Description[gse@result$enrichmentScore<0][1:20]
  gse_pos<-gse@result$Description[gse@result$enrichmentScore>0][1:20]
  
  gse_top_10_each<-c(gse_neg[1:10],gse_pos[1:10])
  
  #setup for emmap plots
  gse_emap<-gse@result[,c('ID','Description','enrichmentScore','pvalue')]
  gse_emap<-gse_emap[1:50,]
  my_gene_set <- prep_gene_set(myGENE2GO,gse_emap,organism_obj)
  termSim_table <- create_go_termsim_table(gse_emap, organism_str, ontology)
  
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/emmap_custom_enhanced.png"),width=50,height=50,
                      plot=emapplot_visEAGO(y_custom=gse_emap,
                                            showCategory =gse_emap$Description,
                                            simTerm_table=termSim_table,
                                            geneSets=my_gene_set,
                                            color="enrichmentScore"))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/emmap_custom_condensed.png"),
                      plot=emapplot_visEAGO(y_custom=gse_emap,
                                            showCategory =gse_emap$Description,
                                            simTerm_table=termSim_table,
                                            geneSets=my_gene_set,
                                            color="enrichmentScore"))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/dotplot_neg.png"),width = 25, height = 20,
                      dotplot_YL_mod_enrich_legend(gse, showCategory = gse_neg, title = exp_name , split=".sign",color="enrichmentScore") + facet_grid(.~.sign))
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/dotplot_pos.png"),width = 25, height = 20,
                      dotplot_YL_mod_enrich_legend(gse, showCategory = gse_pos, title = exp_name , split=".sign",color="enrichmentScore") + facet_grid(.~.sign))
  
  
  my_gg_save_function(save_path_name = paste0(save_path,"/GSEA/dotplot_top_20.png"),width = 25, height = 20,
                      dotplot_YL_mod_enrich_legend(gse, showCategory = gse_top_10_each, title = exp_name , split=".sign",color="enrichmentScore") + facet_grid(.~.sign))
  
}

# GO heat map functions results -------------------

#' function which extracts GO descriptions of genesets from a GSEA object
#'
#' @param x The GSEA object as produced by ClusterProfiler's \code{gseGO} function
#' @param n A vector of GO IDs to be extracted or an integer for top n amount of GO terms
#' 
#' @return None
#'
#' @examples
#'
#' @export
extract_geneSets <- function(x, n) {
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates concept network paths and prepares a list of data for
#' a subsequent heatmap of GO terms.
#' 
#' The function takes in a list containing experiment names which in turn contain 
#' a vector of GO IDs to be extracted/plotted. The function will create a 'custom_GO_selection'
#' in the clusterProfiler_results folder, then create and save the concept network plots for 
#' the specified GO terms.
#' The function will also store the necessary GO information in a list and return the 
#' information. This will serve as a basis for the creation of a GO heatmap.
#'
#' @param main_list A list of vectors where the name represents a experiment (differential
#' gene expression comparison) and the vectors contain the GO IDs of interest.
#' @param path_to_clusterPro The path to the cluster_profiler_results folder which in
#' turn contains the RData objects read by this function
#' 
#' @return heat_map_list A list containing the extracted gene sets from the RData object
#'
#' @examples
#'
#' @export
create_cnetplot_and_heat_map_list <-function(main_list,path_to_clusterPro,plot_cnet=T){
  save_path<-paste0(path_to_clusterPro,'custom_GO_selection')
  dir.create(save_path)
  heat_map_list<-list()
  for (experiment in names(main_list)){
    load(paste0(path_to_clusterPro,experiment,'/cluster_profiler_objects.RData'))
    exp_save_path<-paste0(save_path,'/',experiment)
    dir.create(exp_save_path)
    
    if(plot_cnet==T){
      my_gg_save_function(save_path_name = paste0(exp_save_path,'/all_cnets_enhanced.png'),width=50,height=50,
                          plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = main_list[[experiment]], node_label='all'))
      my_gg_save_function(save_path_name = paste0(exp_save_path,'/all_cnets_condensed.png'),
                          plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = main_list[[experiment]], node_label='category'))
      for (desc in main_list[[experiment]]){
        my_gg_save_function(save_path_name = paste0(exp_save_path,'/',desc,'_condensed.png'),
                            plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = desc, node_label='gene'))
        
        my_gg_save_function(save_path_name = paste0(exp_save_path,'/',desc,'_enhanced.png'),width=50,height=50,
                            plot=cnetplot(gse, categorySize="pvalue", foldChange=gse@geneList, showCategory = desc, node_label='gene'))
      }
    }

    
    heat_map_list[[experiment]]<-extract_geneSets(gse, main_list[[experiment]])
  }
  save(list = c("heat_map_list"), file = paste0(save_path,"/heat_map_list.RData"))
  return(heat_map_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which adds the number of items in a split vector to the names
#' 
#' A split vector is a vector containing the 'splits' of a heatmap, where the vector
#' would be c('control','control','exp','exp','exp'). The 'split' would occur at the
#' control/exp shift. This function adds the number of appearances of each category of splits.
#' In the provided example, the resulting returned factor would be:
#' c('control (2)','control (2)', 'exp (3)', 'exp (3)', 'exp (3)')
#'
#' @param split_vector A vector containing the string used to split a heatmap
#' into separate categories
#' 
#' @return split_vector The updated split_vector converted to a factor.
#'
#' @examples
#'
#' @export
add_length_to_split_vector <-function(split_vector){
  #Add the number of genes in each GO pathway
  longest_val<-unique(split_vector[nchar(split_vector)==max(nchar(unique(split_vector)))])
  numbered_splits<-c()
  for (GO in unique(split_vector)){
    if (table(split_vector)[[GO]]>1){
      num_to_add<-paste0(' (',table(split_vector)[[GO]],')')
    }else{
      num_to_add<-''
    }
    if (GO == longest_val){#Add extra whitespace at the end to force better width of plot
      new_name<-(paste0(split_vector[split_vector==GO][1],num_to_add,'     '))
    }else{
      new_name<-(paste0(split_vector[split_vector==GO][1],num_to_add))
    }
    
    numbered_vect<-rep(new_name,length(split_vector[split_vector==GO]))
    numbered_splits<-c(numbered_splits,numbered_vect)
  }
  
  split_vector <- factor(numbered_splits, levels=unique(numbered_splits))
  return(split_vector)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates a dataframe of genes involved in specific GO terms.
#' The function takes in the object created by \code{create_cnetplot_and_heat_map_list}
#' or read by an Robject.
#' 
#' This list contains the GOs of interest for specific comparisons along with the 
#' names of the genes involved in each GO. The function will use this information
#' to create a large dataframe containing all the relevant genes, for all relevant
#' experiments. It will also create a named vector which contains the name of the 
#' GO pathway repeated as many times as there are genes for that gene set.
#' 
#' This vector will be used to create 'splits' in the subsequent heatmap
#'
#' @param heat_list A list of experiment names followed by GO IDs of interest
#' along with the genes which are used within that GO ID
#' @param sub_all_counts A dataframe containing all the count values for the 
#' genes of interest for all the samples/count files of interest
#' 
#' @return big_GO_df A large dataframe containing the count information for all
#' the genes for all the selected GO IDs
#' @return named_gene_vect A named vector used to 'split' the heatmap per GO term
#' 
#' @examples
#'
#' @export
create_big_GO_df <-function(heat_list,sub_all_counts){
  #Create list of GO pathways with the genes in each pathway
  #Merge duplicate pathways, include all unique genes found
  my_GO_list<-list()
  for (exp in names(heat_list)){
    for (GO_path in names(heat_list[[exp]])){
      #If the GO pathway is seen twice (or more) include all unique genes from each occurrence
      if(GO_path %in% names(my_GO_list)){
        my_GO_list[[GO_path]]<-unique(c(my_GO_list[[GO_path]],heat_list[[exp]][[GO_path]]))
      }else{
        my_GO_list[[GO_path]]<-heat_list[[exp]][[GO_path]]
      }
    }
  }
  
  #Create dataframes of values for each GO
  #Make the gene_names unique for each GO to allow for duplicate genes
  #This is done by appending _GOn to the end of each gene
  #Then load each dataframe into a big dataframe which will serve for the heatmap
  GO_rename_inc=1
  named_gene_vect<-c()
  big_GO_df<-as.data.frame(NULL)
  for (GO_path in names(my_GO_list)){
    sub_GO_df<-sub_all_counts[my_GO_list[[GO_path]],,drop=F]
    rownames(sub_GO_df)=paste0(rownames(sub_GO_df),'_GO',GO_rename_inc)
    # colnames(sub_GO_df)=colnames(subset_all_counts)
    GO_rename_inc=GO_rename_inc+1#Increment for renaming
    
    
    temp_vect=rep(GO_path,nrow(sub_GO_df))
    names(temp_vect)<-c(rownames(sub_GO_df))
    
    named_gene_vect<-c(named_gene_vect,temp_vect)
    
    #Create the dataframe for the heatmap
    if(nrow(big_GO_df)==0){
      big_GO_df<-sub_GO_df
    }else{
      big_GO_df<-rbind(big_GO_df,sub_GO_df)
    }
  }
  
  big_GO_df<-t(big_GO_df)
  
  return(list(big_GO_df,named_gene_vect))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Function which retrieves the log2FoldChange results for a selected set of genes
#' in a selected set of experiments.
#' 
#' The function requires a search list whose names are the experiments to search for
#' log2FoldChange results
#' The function also takes in a vector of genes to extract.
#'
#' @param DE_path The path to the differential expression results
#' @param search_list A vector whose names will be the experiments to extract the L2FCs from
#' @param genes_to_extract A named vector whose values are the genes for which the 
#' L2FC values will be extracted
#' 
#' @return l2fc_df A dataframe where the L2FC value of each gene is given for each experiment
#' @return named_vect_L2FC A named vector containing the experiment names used to find the L2FC values
#' 
#' @examples
#'
#' @export
get_L2FC_results<-function(DE_path,search_list,genes_to_extract){
  #Start by getting the L2FC for each experiment
  genes_of_interest<-unique(unname(genes_to_extract))
  l2fc_df<-as.data.frame(NULL)
  named_vect_L2FC<-c()
  for (exp in names(search_list)){
    DE_df<-read.csv(paste0(DE_path,'/',exp,'/DE_raw_data.csv'))
    DE_df<-DE_df[DE_df$gene_id %in% genes_of_interest,c('gene_id','log2FoldChange')]
    rownames(DE_df)=DE_df$gene_id
    DE_df<-DE_df[genes_of_interest,]#Sort it so it's the same for every file
    colnames(DE_df)=c('gene_id',exp)
    if (nrow(l2fc_df)==0){
      l2fc_df<-DE_df
    }else{
      l2fc_df<-cbind(l2fc_df,DE_df[exp])
    }
    
    temp_vect<-(exp)
    names(temp_vect)<-exp
    named_vect_L2FC<-c(named_vect_L2FC,temp_vect)
  }

  if (length(names(search_list))==1){
    #Only one name, so one column. Use this filter format to prevent the conversion
    #to a dataframe
    l2fc_df<-l2fc_df[2]
  }else{
    l2fc_df<-l2fc_df[,-1]
  }
  
  
  return(list(l2fc_df,named_vect_L2FC))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which handles the plotting of a GO heatmap
#' 
#' The heatmap contains either the count or L2FC values for the genes
#' of chosen GO IDs for all experiments/comparisons selected.
#' 
#' The heatmap is visually split on every GO and experiment in order to better 
#' differentiate the different segments of the heatmap
#' 
#' 
#' @param heat_df The dataframe as produced by \code{create_big_GO_df}
#' @param col_split_vect A factor used to split the columns (GO IDs)
#' @param row_split_vect A factor used to split the rows (experiments)
#' @param save_path The path where the heatmap will be saved
#' @param legend_values_name The name given to the value legend (counts or L2FC)
#' @param save_name The file name that will be given to the heatmap
#' @param legend_subset_name The legend name given to the experiments (subsets, comparisons etc...)
#' @param log_transform Weather or not to log transform the values
#' @param custom_width A width in cm
#' @param custom_height A height in cm
#' 
#' @return None
#' 
#' @examples
#'
#' @export
plot_big_GO_heat <-function(heat_df,col_split_vect,row_split_vect,save_path,legend_values_name='counts',
                            save_name='big_GO_heat',legend_subset_name='Subsets',log_transform=F,
                            custom_width=20,custom_height=10,custom_range=NULL){
  
  if (log_transform==T){
    heat_df<-log10(heat_df+1)
  }
  
  col_splits<-add_length_to_split_vector(col_split_vect)
  row_splits<-add_length_to_split_vector(row_split_vect)
  
  
  #Create a vector of colors of sufficient size for the number of subsets
  #this vector will contain at most 12 different colors
  my_cols<-c('#1f78b4','#33a02c','#fb9a99','#a6cee3','#e31a1c','#b2df8a','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928','#ffff99')
  while(length(col_splits)>length(my_cols)){
    my_cols<-c(my_cols,my_cols)
  }
  
  
  #Set up the top annotation blocks for the subsets
  fill_set_cols=gpar(fill=my_cols[1:length(unique(col_splits))])
  top_annot = HeatmapAnnotation(foo = anno_block(gp = fill_set_cols))
  
  fill_set_rows=gpar(fill=2:(length(unique(row_splits))+1))
  left_annot = rowAnnotation(foo = anno_block(gp = fill_set_rows))
  
  #Create legends
  #Create the legend for the subsets
  lgd_list = list(
    Legend(labels =unique(as.vector.factor(col_splits)), title = "GO pathways", 
           legend_gp = fill_set_cols),
    Legend(labels = unique(as.vector.factor(row_splits)), title = legend_subset_name, 
           legend_gp = fill_set_rows))
  
  
  #find minimum, mid-range, and maximum values to create and label the legend
  if(is.null(custom_range)==F){
    max_val<-4
    min_val<--4
    mid_val<-0
  }else{
    max_val<-ceiling(max(heat_df))
    min_val<-floor(min(heat_df))
    mid_val<-(max_val+min_val)/2
  }
  
  col_fun = colorRamp2(c(min_val, mid_val ,max_val), c("blue","white", "red"))

  
  my_heat<-Heatmap(heat_df,name=legend_values_name,cluster_rows = T,cluster_columns = T,
                   show_column_names = F,show_row_names = F,row_names_side='left',
                   # col=col_fun,
                   row_split = row_splits,column_split=col_splits,
                   border=T,na_col = 'gray',column_title = NULL,row_title = NULL,
                   top_annotation = top_annot,left_annotation=left_annot,
                   cluster_column_slices=F,cluster_row_slices = F,
                   heatmap_legend_param=list(
                     title = legend_values_name, at = c(min_val,mid_val ,max_val),
                     labels = c(min_val,mid_val ,max_val),
                     legend_height = unit(custom_height*0.2, "in")
                   )
                   )
  
  
  save_name_svg<-paste0(save_path,save_name,'.svg')
  svg(save_name_svg,width=custom_width,height=custom_height)
  draw(my_heat, annotation_legend_list = lgd_list)
  dev.off()
  
  save_name_png<-paste0(save_path,save_name,'.png')
  png(save_name_png,width=custom_width*96,height=custom_height*96)
  draw(my_heat, annotation_legend_list = lgd_list)
  dev.off()
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which plots a heatmap for a single GO. Where the heatmap
#' shows the gene values (counts or L2FC) of that GO for each given experiment/subset.
#' 
#' 
#' @param heat_df The dataframe as produced by \code{create_big_GO_df}
#' @param row_split_vect A factor used to split the rows (experiments)
#' @param save_path The path where the heatmap will be saved
#' @param save_name The file name that will be given to the heatmap
#' @param legend_subset_name The legend name given to the experiments (subsets, comparisons etc...)
#' @param log_transform Weather or not to log transform the values
#' @param custom_width A width in cm
#' @param custom_height A height in cm
#' @param sig_gene_list A list containing the significant genes in each experiment
#' 
#' @return None
#' 
#' @examples
#'
#' @export
plot_single_GO_heat <- function(heat_df,row_split_vect,save_path,save_name,legend_values_name='counts',
                                legend_subset_name='Subsets', custom_width=7.5,custom_height=1.5,sig_gene_list=NULL,
                                custom_range=NULL,show_legend=T){
  row_splits<-add_length_to_split_vector(row_split_vect)
  fill_set_rows=gpar(fill=2:(length(unique(row_splits))+1))
  left_annot = rowAnnotation(foo = anno_block(gp = fill_set_rows))
  
  
  bot_annot_names<-colnames(heat_df)
  GO_num_remove<-strsplit(bot_annot_names[1],'_')[[1]][2] #Identify appended val to remove
  bot_annot_names<-unlist(strsplit(bot_annot_names,paste0('_',GO_num_remove))) #Remove appended val from all genes
  
  bot_annot = HeatmapAnnotation(foo = anno_text(bot_annot_names, gp = gpar(fontsize = 10)))
  
  
  
  #Create the legend for the subsets
  lgd_list = list(
    Legend(labels = unique(as.vector.factor(row_splits)), title = legend_subset_name, 
           legend_gp = fill_set_rows))
  
  #find minimum, mid-range, and maximum values to create and label the legend

  if(is.null(custom_range)==F){
    max_val<-4
    min_val<--4
    mid_val<-0
  }else{
    max_val<-ceiling(max(heat_df))
    min_val<-floor(min(heat_df))
    mid_val<-(max_val+min_val)/2
  }

  col_fun = colorRamp2(c(min_val, mid_val ,max_val), c("blue","white", "red"))
  
  my_heat<-Heatmap(heat_df,name=legend_values_name,cluster_rows = T,cluster_columns = T,
                   show_column_names = F,show_row_names = F,row_names_side='left',
                   row_split = row_splits,
                   col=col_fun,
                   border=T,na_col = 'gray',column_title = NULL,row_title = NULL,
                   left_annotation=left_annot, #bottom_annotation = bot_annot,
                   cluster_column_slices=F,cluster_row_slices = F,
                   heatmap_legend_param=list(
                     title = legend_values_name, at = c(min_val,mid_val,max_val), 
                     labels = c(min_val,mid_val,max_val),
                     legend_height = unit(custom_height*0.4, "in")
                   ),
                   show_heatmap_legend = show_legend,
                   #Function that marks significant genes in the heatmap
                   layer_fun = function(j, i, x, y, w, h, fill) {
                     #Create vector of gene names in clustered order
                     clustered_gene_order<-names(heat_df[i[1],j])
                     clustered_gene_order<-gsub(clustered_gene_order,pattern='_.*',replacement = '')
                     #Find comparison being used/plotted
                     current_comp<-row.names(heat_df)[unique(i)]
                     #Extract significant genes, find corresponding indexes in clustered gene names
                     sig_genes_interest<-sig_gene_list[[current_comp]][[save_name]]
                     gene_marker_index=which(clustered_gene_order %in% sig_genes_interest)
                     #Create index matrix for grid_points
                     ind_mat = restore_matrix(j, i, x, y)
                     #If there are indexes for significant genes, add them to heatmap
                     if (length(gene_marker_index)>0){
                       grid.points(x[gene_marker_index], y[gene_marker_index], pch = 16, size = unit(2, "mm"))
                     }
                     
                   })
  
  #Extends the width of the svg based on the number of genes
  #This enables the svg to be used as a means to identify location of specific genes
  extend_val_for_svg=round(length(bot_annot_names)/100)
  if (extend_val_for_svg==0){
    svg_width=custom_width
  }else{
    svg_width=custom_width*extend_val_for_svg
  }
  
  if(show_legend==T){
    save_name_svg<-paste0(save_path,save_name,'.svg')
    svg(save_name_svg,width=svg_width,height=custom_height)
    draw(my_heat, annotation_legend_list = lgd_list,column_title =save_name)
    dev.off()
    
    save_name_png<-paste0(save_path,save_name,'.png')
    png(save_name_png,width=custom_width*96,height=custom_height*96)
    draw(my_heat, annotation_legend_list = lgd_list,column_title =save_name)
    dev.off()
  }else{
    save_name_svg<-paste0(save_path,save_name,'.svg')
    svg(save_name_svg,width=svg_width,height=custom_height)
    draw(my_heat,column_title =save_name)
    dev.off()
    
    save_name_png<-paste0(save_path,save_name,'.png')
    png(save_name_png,width=custom_width*96,height=custom_height*96)
    draw(my_heat,column_title =save_name)
    dev.off()
  }
  

  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Wrapper function which will create a single heatmap per GO ID. Each heatmap will contain the gene
#' information for all the provided subsets/expeirments.
#' 
#' 
#' 
#' @param heat_df The dataframe as produced by \code{create_big_GO_df}
#' @param GO_gene_vect Named vector where names are GO IDs and values are gene names
#' @param row_name_vector Named vector where names are samples/count files and 
#' values are experiments/subsets
#' @param save_path The path where the heatmaps will be saved
#' @param folder_name The folder that will be created to store the heatmaps
#' @param use_L2FC Boolean indicating if L2FCs are being plotted, used for legend title
#' @param sig_gene_list List of comparisons and GOs containing significant genes
#' 
#' @return None
#' 
#' @examples
#'
#' @export

plot_single_heats_wrapper <- function(heat_df,GO_gene_vect,row_name_vector,save_path,folder_name,use_L2FC=F,sig_gene_list=NULL,
                                      custom_range=NULL,custom_width=7.5,custom_height=1.5,show_legend=T){
  #Create folder for single GO heatmaps
  save_path<-paste0(save_path,'custom_GO_selection')
  small_heat_save_path<-paste0(save_path,'/',folder_name)
  
  #Adjust the legend and NULL the significant gene list if counts are being used
  #In a count heatmap, there can be no 'significant' genes as individual samples are 
  #plotted, not comparisons
  if(use_L2FC==T){
    legend_text='FC'
    sig_gene_list=sig_gene_list
  }else{
    legend_text='Counts'
    sig_gene_list=NULL
  }
  
  #Only create the plots if there is nothing in the folder or if the folder does not exist.
  if (!small_heat_save_path %in% list.dirs(save_path)){
    small_heat_save_path<-paste0(small_heat_save_path,'/')
    dir.create(small_heat_save_path)
    #Iterate over each GO and create single GO heatmaps
    for (i in 1:length(unique(unname(GO_gene_vect)))){
      sub_gene_vect<-GO_gene_vect[endsWith(names(GO_gene_vect),paste0('_GO',i))]
      sub_bg_GO<-heat_df[,names(sub_gene_vect)]
      GO_name<-unique(sub_gene_vect)
      
      plot_single_GO_heat(sub_bg_GO,row_name_vector,save_path = small_heat_save_path, 
                          save_name=GO_name,legend_values_name = legend_text,
                          sig_gene_list=sig_gene_list,
                          custom_range=custom_range,
                          custom_width = custom_width,custom_height = custom_height,
                          show_legend=show_legend)
      
    }
  }else{
    message('The inputed folder name already exists in the provided location')
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Function that creates a excel file for each comparison
#' 
#' The excel file contains the differential gene expression results for the genes
#' of the GOs in the provided list (exp_GO_list). Each GO is given a separate sheet
#' of the excel. An additional column is added where a YES or NO is given in regards
#' to that genes significance using the log2foldchange and either the pvalue or adjusted pvalue.
#' Respectively these filters are greater/lower than +1/-1 and smaller than 0.05.
#' 
#' The excel files are saved to the indicated location and they carry the name of the
#' comparison they represent.
#' 
#' @param exp_GO_list A nested list containing comparisons, GOs, and genes (in that order)
#' @param path_DE The path to the differential gene expression files
#' @param significance_p The pvalue type to use for labelling significant genes
#' can be either pvalue or padj (adjusted p value)
#' @param excel_save_path The location where the excel files will be saved
#' 
#' @return significant_genes A list of significant gene_IDs for each comparison
#' 
#' @examples
#'
#' @export
create_excel_data_files<-function(exp_GO_list,path_DE,significance_p,excel_save_path){
  dir.create(excel_save_path)
  significance_p='padj'
  sig_col_name=paste0('significant_',significance_p)
  library(readxl)
  significant_genes<-list()
  for (exp in names(exp_GO_list)){
    print(exp)
    significant_genes[[exp]]<-list()
    excel_name<-paste0(exp,'.xlsx')
    excel_file<-paste0(excel_save_path,excel_name)
    for (GO in names(exp_GO_list[[exp]])){
      print('##########')
      print(GO)
      print('##########')
      genes_of_GO<-exp_GO_list[[exp]][[GO]]
      DE_file<-read_csv(paste0(path_DE,'/',exp,'/DE_raw_data.csv'))
      DE_file<-DE_file[DE_file$gene_id %in% genes_of_GO,]
      DE_file[[sig_col_name]] <- "NO"
      DE_file[[sig_col_name]][DE_file[[significance_p]] <0.05 & abs(DE_file$log2FoldChange)>1] <- "YES"
      sig_genes<-DE_file$gene_id[DE_file[[sig_col_name]]=="YES"]
      DE_file<-as.data.frame(DE_file)
      significant_genes[[exp]][[GO]]<-sig_genes
      if (excel_name %in% list.files(excel_save_path)==F){
        file <- paste(excel_file, sep = "")
        write.xlsx(DE_file, file, sheetName = GO,row.names = F) 
      }else{
        sheet_names<-excel_sheets(path = excel_file)
        sub_GO_name<-paste(strsplit(GO,'')[[1]][1:31],collapse='')
        if(sub_GO_name %in% sheet_names==T){
          sheet_name<-strsplit(GO,'')[[1]][1:30]
          sheet_name<-c(sheet_name,1)
          sheet_name<-paste(sheet_name, collapse = "")
        }else{
          sheet_name<-GO
        }
        write.xlsx(DE_file, file, sheetName = sheet_name, append = TRUE,row.names = F) 
      }
    }
    
  }
  return (significant_genes)
}


# Process visEAGO results -------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A wrapper function which processes visEAGO results in a format which can be 
#' used in clusterProfiler plotting.
#' 
#' @param vis_clust_data The visEAGO object as produced by the visEAGO pipeline
#' @param cluster_number The cluster of interest
#' @param target_exp The experiment of interest
#' @param organism_str A string of the annotation dbi organism to be used
#' @param organism_obj The annotation dbi object ot be used.
#' @param ontology The ontology to be used
#' 
#' @return my_termSim The table of term similarities
#' 
#' @example 
#' 
#' 
#' 
#' @export
process_visEAGO_results_for_visu <- function(vis_clust_data,cluster_number,target_exp,organism_str,organism_obj,ontology='BP'){
  
  #Set up the GENE2GO annotation
  Bioconductor<-ViSEAGO::Bioconductor2GO()
  myGENE2GO<-ViSEAGO::annotate(
    organism_str,
    Bioconductor)
  G2GO <- myGENE2GO@BP
  
  #Filter and format the visEAGO object to obtain the GOs from the cluster
  #GOs which are also found in the GENE2GO object are stored in 'clus_df_complete'
  #While GOs which were in the cluster but not in GENE2GO are stored in 'GO_not_found'
  not_found_found_GOs <- prep_format_clustered_GOs(vis_clust_data,G2GO,target_exp,organism_obj=organism_obj,cluster_number=cluster_number)
  
  go_not_found <- not_found_found_GOs[[1]]
  clus_df_complete <- not_found_found_GOs[[2]]
  
  #Check and update go_found/cluster_df_complete
  #Sometimes genes are duplicated and therefore the total gene count of the GOs are off
  #This function checks for such events and updates the object
  clus_df_complete <- check_update_go_found(clus_df_complete)
  
  #Creates a geneset and termSim table using all the GOs found
  my_gene_set <- prep_gene_set(myGENE2GO,clus_df_complete,organism_obj)
  termSim_table <- create_go_termsim_table(clus_df_complete, organism_str, ontology)
  
  #Prepares a list with all the objects to be returned.
  return_list <- list(go_not_found,clus_df_complete,my_gene_set,termSim_table)
  
  return (return_list)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A function which converts a visEAGO cluster object to a data object that can be 
#' used to plot the pathways found. 
#' 
#' visEAGO splits it's results into clusters, this function retrieves the GO pathways
#' in this cluster, it then attempts to complete the object by finding the genes in 
#' each pathway using the G2GO object (bioconducter GO annotation). GO pathways which were found
#' in the G2GO pathway are processed and stored in 'clus_df_complete', while GO pathways
#' which were in the visEAGO cluster but not in the bioconductor annotation are stored
#' in 'go_not_found', these will be processed in a later part of the pipeline.
#' 
#' @param vis_clust_data The visEAGO object to be processed
#' @param G2GO A dataframe version of the GENE2GO_object for a specific ontology
#' @param exp_of_interest The experiment of interest
#' @param organism_obj A annotation dbi annotation object for the target organism
#' @param cluster_number The cluster of interest
#' 
#' @return return_list
#' 
#' @example 
#' 
#' 
#' 
#' @export
prep_format_clustered_GOs <- function(vis_clust_data,G2GO,exp_of_interest,organism_obj,cluster_number){
  #Convert the visEAGO enriched go data to a dataframe
  full_df <- as.data.frame(vis_clust_data@enrich_GOs@data)
  
  #Columns 1 to 5 are as follows:
  #GO.cluster, GO.ID, term, definition, pvalue
  cols_to_keep <- (colnames(full_df)[1:5])
  #The visEAGO object contains the results for several experiments, here we isolate the 
  #result for the specific experiment of interest (we remove columns which do not correspond to
  #the experiment of interest)
  cols_to_keep <- c(cols_to_keep,colnames(full_df)[startsWith(colnames(full_df), exp_of_interest)])
  
  #Perform the filter for both the experiment and the cluster number of interest
  exp_df <- full_df[cols_to_keep]
  exp_clus <- exp_df[exp_df$GO.cluster==cluster_number,]
  
  
  #Some column renaming and data formating
  exp_clus_nams <- c()
  for (nam in colnames(exp_clus)){
    if (startsWith(nam, exp_of_interest)==T){
      exp_clus_nams <- c(exp_clus_nams, str_remove(nam,paste0(exp_of_interest,'.')))
    }else{
      exp_clus_nams <- c(exp_clus_nams,nam)
    }
  }
  colnames(exp_clus)=exp_clus_nams
  mod_clus <- exp_clus %>% 
    dplyr::select(c('GO.ID','term','genes_frequency','pvalue','Significant_genes_symbol')) %>% 
    dplyr::rename(ID = GO.ID, Description=term, GeneRatio=genes_frequency)
  mod_clus$Significant_genes_symbol <- gsub(';', '/', mod_clus$Significant_genes_symbol)
  
  #Below we create a dataframe of the G2GO annotation for the GOs present in the
  #visEAGO cluster of interest. Genes in G2GO are stored as 
  #EntrezID, therefore a conversion is done to convert them to symbol.
  entrez_vect <- c()
  GO_vect <- c()
  for (entrez in names(G2GO)){
    for (GO_id in (unname(G2GO[[entrez]]))){
      if (GO_id %in% mod_clus$ID){
        entrez_vect <- c(entrez_vect,entrez)
        GO_vect <- c(GO_vect,GO_id)
      }
    }
  }
  gene_symb_vect <- as.vector(mapIds(organism_obj, keys = entrez_vect, keytype = "ENTREZID", column="SYMBOL"))
  GO_in_data <- data.frame(ID=GO_vect, GENE_SYMBOL=gene_symb_vect)
  
  
  #Some formatting and creation of the dataframe
  geneID_vector <- c()
  count_vector <- c()
  for (GO_ID in GO_in_data$ID){
    sub_GO <-GO_in_data[GO_in_data$ID==GO_ID,]
    found_count <- length(sub_GO$GENE_SYMBOL)
    geneID_vector <- c(geneID_vector,paste(sub_GO$GENE_SYMBOL,collapse="/"))
    count_vector <- c(count_vector,found_count)
  }
  GO_in_data$geneID <- geneID_vector
  GO_in_data$Count <- count_vector
  GO_in_data <- GO_in_data[c('ID','geneID','Count')]
  GO_in_data <- unique(GO_in_data)
  
  
  #Here we split by the GOs which were found in G2GO and those who were not
  mod_clus = subset(mod_clus, select = -c(GeneRatio) )
  go_not_found <- mod_clus[!(mod_clus$ID %in% GO_in_data$ID),]
  
  #For the GOs found, we merge it with the annotation data
  clus_df_complete <- merge(mod_clus,GO_in_data, by='ID')
  rownames(clus_df_complete)=clus_df_complete$ID
  
  #Store results in list format to be able to return both
  return_list <- list(go_not_found,clus_df_complete)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which updates the go_found object
#' 
#' First checks that the genes which are significant are also counted as genes
#' Which are part of the pathway, NAs are removed, and duplicates are also removed.
#' 
#' If nay of the above occurred, the 'count' column is updated as the number of
#' genes in the pathway may have changed as a result of the update.
#' 
#' @param go_found a dataframe of found GOs
#' 
#' @return my_gene_set A dataframe for the gene set
#' 
#' @example 
#' 
#' 
#' 
#' @export
check_update_go_found <- function(go_found){
  for (i in nrow(go_found)){
    update=F
    sig_gene_vect <- as.vector(strsplit(go_found$Significant_genes_symbol[i], '/')[[1]])
    genes_in_go <- as.vector(strsplit(go_found$geneID[i], '/')[[1]])
    
    
    for (gene in sig_gene_vect){
      if (!gene %in% genes_in_go){
        genes_in_go <- c(genes_in_go,gene)
        update=T
      }
    }
    genes_in_go<-genes_in_go[genes_in_go!="NA"]
    genes_in_go<-genes_in_go[!is.na(genes_in_go)]
    genes_in_go<-unique(genes_in_go)
    if (update==T){
      go_found$Count[i] <- length(genes_in_go)
      go_found$geneID[i] <- paste(genes_in_go,collapse="/")
    }
  }
  
  return(go_found)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which creates a gene set for a set of custom GO pathways.
#' 
#' A gene set is a dataframe which provides the genes in each GO pathway provided
#' The function utilizes a GENE2GO object to create a gene set for the custom
#' GO pathways. The function is specific to BP, but is easily swapped.
#' 
#' The genes are initially found as EntrezID but are converted to Symbol ID
#' 
#' @param GENE2GO_object The GENE2GO object from which the gene set can be extracted
#' @param go_found A dataframe containing a 'ID' column with the GO pathways to be used
#' @param organism The GO.DB object for the organism
#' 
#' @return my_gene_set A dataframe for the gene set
#' 
#' @example 
#' 
#' 
#' 
#' @export
prep_gene_set <- function(GENE2GO_object,go_found,organism){
  G2GO <- GENE2GO_object@BP
  my_gene_set <- list()
  for (entrez in names(G2GO)){
    for (go_id in G2GO[[entrez]]){
      if (go_id %in% go_found$ID){
        my_gene_set[[go_id]]<-c(my_gene_set[[go_id]],entrez)
      }
    }
  }
  
  
  for (go_id in names(my_gene_set)){
    entrez_vect <-my_gene_set[[go_id]] 
    entrez_vect <- suppressMessages(as.vector(mapIds(organism, keys = entrez_vect, keytype = "ENTREZID", column="SYMBOL")))
    my_gene_set[[go_id]]<-unique(entrez_vect)
  }
  
  return (my_gene_set)
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which creates the terminology similarity (termsim) table for 
#' a set of custom GO pathways.
#' 
#' Function uses a GOSEMSIM innate function to create the term sim table
#' for the GOs provided. The GO terms should be in a dataframe within an 
#' 'ID' column
#' 
#' @param found_GOs A dataframe containing a column ID with the GO terms
#' @param organism_str The GO.DB organism in string format
#' @param ontology The ontology to be used
#' 
#' @return my_termSim The table of term similarities
#' 
#' @example 
#' 
#' 
#' 
#' @export
create_go_termsim_table <- function(found_GOs, organism_str,ontology){
  go_IDs <- found_GOs$ID
  
  #Create go_data for organism
  d <- godata(organism_str, ont=ontology)
  #Create term similarity table
  my_termSim <- round(termSim(go_IDs, go_IDs, d),digits=3)
  my_termSim <- melt(my_termSim)
  my_termSim <- my_termSim[my_termSim[,1] != my_termSim[,2],]
  my_termSim <- my_termSim[!is.na(my_termSim[,3]),]
  
  return(my_termSim)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which first checks if there are any GO's which were not found
#' in the G2GO annotation, this is translated by the precense of GOs in the 
#' 'GO_not_found' object. If all GO's were found (0 rows) then the function
#' states that all GOs were found and returns the inputed objects as they were.
#' 
#' If some GOs were not found, it will first create a 'missing_GOs' folder in the
#' save path, and search this folder for a file which will supplement for the missing GOs
#' 
#' #If the necessary file is not there, a message is printed to instruct the user with what
#' needs to be done.
#' 
#' If the necessary file is there, a function is called to update the 'GO_found' object with
#' the necessary information to add the GOs which were previously not found.
#' 
#' @param vis_cluster_res_list a list of processed visEAGO objects as produced 
#' by process_visEAGO_results_for_visu()
#' @param target_path path to which the missing GOs folder will be created/searched
#' @param organism_str the string of the annotation db organism
#' @param ontology the ontology (often BP)
#' 
#' @return vis_cluster_res_list
#' 
#' @example 
#' 
#' 
#' 
#' @export
check_missing_GOs <- function(vis_cluster_res_list,target_path,organism_str, ontology){
  
  go_not_found <- vis_cluster_res_list[[1]]
  go_found <- vis_cluster_res_list[[2]]
  my_gene_set <- vis_cluster_res_list[[3]]
  my_termSim <- vis_cluster_res_list[[4]]
  
  missing_desc <- c()
  bool_print_message=F
  if (nrow(go_not_found)>0){
    if (!'missing_GOs' %in% list.files(target_path)){
      bool_print_message =T
    }else{
      for (desc in go_not_found$Description){
        if (! paste0(desc,'.tsv') %in% list.files(paste0(target_path,'/missing_GOs'))){
          missing_desc <-c(missing_desc,desc)
          bool_print_message =T
        }
      }
    }
  
    if (bool_print_message==T){
      dir.create(paste0(target_path,'/missing_GOs'))
      print("There are missing GOs form the analysis, to supplement the analysis, do the following:")
      print(paste0("Use the QuickGo website to search for the missing GO term (",go_not_found$ID,")"))
      print("Make sure to use the proper filters (organisms), then extract the results as a tsv using the export button")
      print("Do the above for each missing go, rename the extracted file to the description of the missing GOs")
      print("put the files inthe 'missing_GOs' folder that was created in the save path of this script")
      print(paste0("The descriptions are: ",go_not_found$Description))
      return()
    }else{
      #Getting here means we have what we need, so we call a function that will
      #Perform the necessary adjustments to the GO_not_found and add it to go_found
      new_results <- fix_results_with_quickGO(vis_cluster_res_list,target_path, organism_str, ontology)
      
      return(new_results)
    }
  }else{
    print('All GOs were found')
    return(vis_cluster_res_list)
  }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' A function which updates our found GOs with information from QuickGO
#' 
#' In the event that a GO pathway was not found in the GENE2GO annotation, a quickGO
#' file can be downloaded by the user and added to the 'missing_GOs' folder. That file
#' is the used in this function to update the missing GOs with the information within the file.
#' 
#' The function first reads the file, filters it for UniProtKB gene product,
#' and removes any duplicate gene symbol. Some data formatting is then performed
#' so that the updated go_found object can be used for plotting.
#' 
#' A new geneset and termsim table is created with the now updated go_found object
#' 
#' @param vis_cluster_res_list a list of processed visEAGO objects as produced 
#' by process_visEAGO_results_for_visu()
#' @param target_path path to which the missing GOs folder will be created/searched
#' @param organism_str the string of the annotation db organism
#' @param ontology the ontology (often BP)
#' 
#' @return return_list A list of four objects - go_not_found, go_found, my_gene_set, my_termsim
#' 
#' @example 
#' 
#' 
#' 
#' @export
fix_results_with_quickGO <- function(vis_cluster_res_list,target_path, organism_str, ontology){
  go_not_found <- vis_cluster_res_list[[1]]
  go_found <- vis_cluster_res_list[[2]]
  my_gene_set <- vis_cluster_res_list[[3]]
  my_termSim <- vis_cluster_res_list[[4]]
  
  
  for (i in 1:nrow(go_not_found)){
    quickGO_res <- read_tsv(paste0(target_path,'/missing_GOs/',go_not_found$Description[i],'.tsv'))
    quickGO_res <- quickGO_res[quickGO_res$`GENE PRODUCT DB`=='UniProtKB',]
    unique_genes <- unique(quickGO_res$SYMBOL)
    updated_not_found <- go_not_found[i,]
    updated_not_found$geneID <- paste(unique_genes,collapse="/")
    updated_not_found$Count <- length(unique_genes)
    rownames(updated_not_found)=updated_not_found$ID
    go_found <- rbind(go_found,updated_not_found)
    
    #Double check that go_found is correct, update if not
    go_found <- check_update_go_found(go_found)
    
    my_gene_set[[go_not_found$ID[i]]] <- unique_genes
  }
  
  
  
  my_termSim <- create_go_termsim_table(go_found, organism_str, ontology)
  
  return_list <- list(go_not_found,go_found,my_gene_set,my_termSim)
  
  return(return_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which creates a single csv per GO pathway.
#' 
#' For each GO pathway, the function retrieves the genes within the pathway, associated
#' the differential gene expression data for those genes (if there is data, otherwise NA),
#' it also provides a description of the gene. The description of the gene in regards to the 
#' pathway is obtained from a gaf (gene annotation file). Lastly the file indicates if the 
#' genes are significant in regards to their differential gene expression results.
#' 
#' 
#' 
#' @param go_found Object containing the GOs of interest
#' @param DE_file The differential gene expresison results for the target experiment
#' @param gaf_file The gaf file
#' @param cluster_save_path The path in which to save the resulting CSVs
#' @param p_filter The filter to use for significance, expects adjustedPvalue or pvalue
#' 
#' @return None
#' 
#' @example 
#' 
#' 
#' 
#' @export
create_go_gene_files <- function(go_found,DE_file,gaf_file,cluster_save_path,p_filter='pvalue'){
  dir.create(paste0(cluster_save_path,'/go_gene_files'))
  go_save_path <- paste0(cluster_save_path,'/go_gene_files')
  for (i in 1:nrow(go_found)){
    all_genes_in_path <- strsplit(go_found$geneID[i], '/')[[1]]
    # all_genes_in_path <- all_genes_in_path[all_genes_in_path!="NA"]
    
    
    DE_genes <- DE_file[DE_file$gene_id %in% all_genes_in_path ,]
    DE_genes <- DE_genes[,1:7]
    
    # go_annot <- read_tsv('mgi.gaf', skip = 52, col_names = FALSE, guess_max=2000)
    
    go_annot <- gaf_file
    
    go_annot <- go_annot[go_annot$X3 %in% all_genes_in_path,]
    go_annot <- go_annot[!duplicated(go_annot$X3),]
    go_annot <- go_annot[,c('X3','X10')]
    colnames(go_annot)=c('gene_id','Description')
    
    go_file <- as.data.frame(all_genes_in_path)
    colnames(go_file)=c('gene_id')
    
    go_file <- merge(go_file,go_annot,by='gene_id',all = TRUE)
    go_file <- merge(go_file,DE_genes,by='gene_id',all = TRUE)
    
    sig_vector<-c()
    for (ii in 1:nrow(go_file)){
      if (is.na(go_file$log2FoldChange[ii])==F){
        if (abs(go_file$log2FoldChange[ii])>1 && go_file[p_filter][ii,] < 0.05){
          sig_vector <- c(sig_vector,'yes')
          next
        }
      }
      sig_vector <- c(sig_vector,'no')
    }
    go_file <- add_column(go_file, significance = sig_vector, .after = "Description")
    
    write.csv(go_file,file = paste0(go_save_path,'/',go_found$Description[i],'.csv'),row.names = F)
  }
}


#' A function which splitsa go term into several possible types. Then plots a cnetplot
#' 
#' Depending on the type selected, the function filters for the selected gene type
#' which is either down (<0), down_sig (<=-1), up (>0), or up_sig (>1). It then plots
#' the resulting value in a newly created folder within the intended location.
#' 
#' 
#' @param gene_set_sub a list with the description as the key and genes within the pathway
#' as the values.
#' @param log2FC_gene_list A list of genes with associated log2FoldChange values
#' @param the_go_id the object which contains the description of the GO term
#' @param clust_save_path The path in which to save the results
#' @param type The type of subsetting to be done
#' 
#' @return None
#' 
#' @example 
#' 
#' 
#' 
#' @export
split_GO_upreg_downreg <- function(gene_set_sub,log2FC_gene_list,the_go_id,clust_save_path,type){
  
  if (type=='down'){
    new_gene_list <-log2FC_gene_list[log2FC_gene_list<0]
  }else if(type=='down_sig'){
    new_gene_list <-log2FC_gene_list[log2FC_gene_list<=-1]
  }else if(type=='up'){
    new_gene_list <-log2FC_gene_list[log2FC_gene_list>0]
  }else if(type=='up_sig'){
    new_gene_list <-log2FC_gene_list[log2FC_gene_list>=1]
  }else{
    print('Not a valid type input')
  }
  
  new_path=paste0(clust_save_path,'/cnetplots/split_plots/')
  dir.create(new_path, showWarnings = FALSE)
  
  name_term=names(gene_set_sub)
  new_gene_sub<-list()
  for (go_name in name_term){
    temp_hold<-list()
    temp_hold[[go_name]]<- gene_set_sub[[go_name]][gene_set_sub[[go_name]] %in% names(new_gene_list)]
    
    my_gg_save_function(paste0(new_path,names(temp_hold),'_',type,'_condensed.png'),
                        plot=cnetplot_mod(temp_hold, categorySize="pvalue",
                                          foldChange=new_gene_list,
                                          showCategory = the_go_id$Description,
                                          node_label="gene"))
    
    my_gg_save_function(paste0(new_path,names(temp_hold),'_',type,'_enhanced.png'),
                        plot=cnetplot_mod(temp_hold, categorySize="pvalue",
                                          foldChange=new_gene_list,
                                          showCategory = the_go_id$Description,
                                          node_label="gene"),width=50,height=50)
    
    
    new_gene_sub[[go_name]]<-temp_hold[[go_name]]
  }
  
  if (length(names(new_gene_sub))>1){
    my_gg_save_function(paste0(new_path,'merged_',type,'_condensed.png'),
                        plot=cnetplot_mod(new_gene_sub, categorySize="pvalue",
                                          foldChange=new_gene_list,
                                          showCategory = the_go_id$Description,
                                          node_label="gene"))
    
    my_gg_save_function(paste0(new_path,'merged_',type,'_enhanced.png'),
                        plot=cnetplot_mod(new_gene_sub, categorySize="pvalue",
                                          foldChange=new_gene_list,
                                          showCategory = the_go_id$Description,
                                          node_label="gene"),width=50,height=50)
  }

  
}


# Modified ClusterProfiler functions plotting functions -------------------

# For the purpose of custom plots, certain functions had to be slightly adjusted, 
# while others needed only be present in this utility script.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A modified version of 'get_igraph' from the enrichplot package 
#' 
#' @param simTerm_table The termsim table
#' @param geneSets The geneSets
#' @param y Dataframe containing the descriptions of the GOs to be plotted
#' @param n Number of enriched terms to display.
#' @param color color variable that used to color enriched terms, e.g. 'pvalue',
#' 'p.adjust' or 'qvalue'.
#' @param cex_line Scale of line width.
#' @param min_edgeThe minimum similarity threshold for whether 
#' two nodes are connected, should between 0 and 1, default value is 0.2.
#' 
#' 
#' @return 
#' 
#' @example 
#' 
#' 
#' 
#' @export
get_igraph_mod <- function(simTerm_table,geneSets,y,  n, color, cex_line, min_edge){
  
  y <- y[match(n, y$Description),]
  n <- length(n)
  g <- emap_graph_build_mod(y = y, simTerm_table = simTerm_table,geneSets=geneSets, color = color,
                            cex_line = cex_line, min_edge = min_edge,
                            pair_sim = NULL, method = NULL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A modified version of 'emap_graph_build' from enrichplot package
#' 
#' @param y a data.frame of GO results
#' @param simTerm_table A semantic similarity table
#' @param geneSets a list gene sets with the names of enrichment IDs
#' @param color a string, the column name of y for nodes colors
#' @param cex_line scale of line width
#' @param min_edge minimum percentage of overlap genes to display the edge,
#' should between 0 and 1, default value is 0.2
#' @param method method of calculating the similarity between nodes,
#' 
#' @return 
#' 
#' @example 
#' 
#' 
#' 
#' @export
emap_graph_build_mod <- function(y, simTerm_table,geneSets, color, cex_line, min_edge,
                                 pair_sim = NULL, method = NULL) {
  wd<- simTerm_table
  if (is.null(method) ==T) {
    # map id to names
    wd[, 1] <- y[wd[, 1], "Description"]
    wd[, 2] <- y[wd[, 2], "Description"]
  }
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
  cnt <- sapply(geneSets[idx], length)
  V(g)$size <- cnt
  colVar <- y[idx, color]
  # V(g)$name <- as.character(y$Description)
  V(g)$color <- colVar
  return(g)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A replicate of the 'add_node_label' of the enrichplot package
#' 
#' The function was added here to allow modified functions access to this function
#' 
#' @param p a ggplot2 object.
#' @param data it is uesd as the `data` parameter of function `ggraph::geom_node_text`, a data.frame or NULL.
#' @param label_location a data.frame with the location of group label.
#' @param label_size_node a numeric value to indicate the font size of the node label.
#' @param cex_label_node a numeric value to indicate the scale of node label size.
#' @param shadowtext  a logical value, whether to use shadow font. 
#' @return a ggplot2 object.
#' 
#' @example 
#' 
#' 
#' 
#' @export
add_node_label_mod <- function(p, data, label_size_node, cex_label_node, shadowtext) {
  if (shadowtext) {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            size = label_size_node * cex_label_node, bg.color = "white", repel=TRUE)
  } else {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            size = label_size_node * cex_label_node, repel=TRUE)
  }
  return(p)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Convert a list of gene IDs to data.frame object.
#' Function originates from the enrichplot package, it was recreated here
#' to allow modified functions to access it
#'
#'
#' @title Convert gene IDs to data.frame object
#' @param inputList A list of gene IDs
#' 
#' @return a data.frame object.
#' 
#' @example 
#' 
#' 
#' 
#' @export
list2df_mod <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Convert a list of gene IDs to igraph object.
#' Function originates from the enrichplot package, it was recreated here
#' to allow modified functions to access it
#'
#' @title Convert gene IDs to igraph object
#' @param inputList A list of gene IDs.
#' @return A igraph object.
#' 
#' @example 
#' 
#' 
#' 
#' @export
list2graph_mod <- function(inputList) {
  x <- list2df_mod(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}

