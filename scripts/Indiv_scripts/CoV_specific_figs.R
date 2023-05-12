setwd('~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/')
source('R/my_WGCNA_functions.R')

load('WGCNA_results/WGCNA_part2.RData')

colnames(datTraits)=gsub(x = colnames(datTraits),pattern = '_',replacement = ' ')

colnames(corr_mat)=gsub(x = colnames(corr_mat),pattern = '_',replacement = ' ')
row.names(corr_mat)=gsub(x = row.names(corr_mat),pattern = '_',replacement = ' ')

library(corrplot)
plot_clinical_corr_plot_V2(datTraits,corr_mat,my_format_p)

svg(paste0("WGCNA_results/clinical_values_heatmap.svg"),15,10)
print(plot_clinical_heatmap(datTraits,traitColors,'Critical'))
dev.off()


#DOTPLOTS!!!!

library(scales)
dotplot_YL_custom <- function(object, x = "geneRatio", color = "p.adjust",
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


source('R/my_ClusterProfiler_functions.R')

path<-'cluster_profiler_results/'
list.files(path)

for(exp in list.files(path)){
  load(paste0(path,exp,'/cluster_profiler_objects.RData'))
  rm(go_enrich_down,go_enrich_up)
  
  desc_list<-c(gse@result$Description[gse@result$enrichmentScore>0][1:10],gse@result$Description[gse@result$enrichmentScore<0][1:10])
  
  my_gg_save_function(save_path_name = paste0(exp,'.png'),width = 26, height = 16,
                      dotplot_YL_custom(gse, label_format = 60, showCategory = desc_list, title = exp , split=".sign",color="enrichmentScore") + facet_grid(.~.sign))
  
  }


