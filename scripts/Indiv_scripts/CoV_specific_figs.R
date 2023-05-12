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


plot_cluster_traj_custom<-function(object,ts_data,ts_mean_data,num_col=4,rem_legend_axis=F){
  if('log10_timepoint' %in% colnames(ts_data)){
    plt <- ggplot(ts_data, aes(y = trans_mean , x = `log10_timepoint`, color = group))
  }else{
    plt <- ggplot(ts_data, aes(y = trans_mean , x = timepoint, color = group))
  }
  plt <- plt + scale_color_manual(values=object@group_colors) +
    geom_line(aes(group = gene_id), alpha = 0.4) +
    geom_point() +
    geom_line(
      data = ts_mean_data, lwd = 1.5, color = "grey50",
      aes(group = group)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0))+
    facet_wrap(~labels, scales = 'free_x', ncol = num_col)+
    theme(strip.text.x = element_text(size = 13))
  
  if(rem_legend_axis==T){
    plt<-plt + theme(legend.position = "none") +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank())
      
  }
  
  return(plt)
}


wrapper_cluster_trajectory_custom<-function(object,cluster_traj_dta,mean_cluster_traj_dta,log_TP=F,plot_name='Ctraj'){
  clust_order<-unique(cluster_traj_dta[,c('cluster','nGenes')])
  clust_order<-clust_order$cluster[order(-clust_order$nGenes)]
  num_needed_figures<-ceiling(length(clust_order)/8)
  #Iterate over number of necessary figures
  for (idx in 1:num_needed_figures){
    #Filter the cluster_map dataframe for the required clusters
    max_clust=8*idx
    if (idx==1){
      min_clust<-1
    }else{
      min_clust<-8*(idx-1)+1
    }
    clusters_to_plot<-clust_order[min_clust:max_clust]
    clusters_to_plot<-clusters_to_plot[!is.na(clusters_to_plot)]#Remove NAs
    
    sub_ts_data<-cluster_traj_dta[cluster_traj_dta$cluster %in% clusters_to_plot,]
    sub_ts_data<-sub_ts_data[order(match(sub_ts_data$cluster,clusters_to_plot)),]
    
    sub_ts_means<-mean_cluster_traj_dta[mean_cluster_traj_dta$cluster %in% clusters_to_plot,]
    sub_ts_means<-sub_ts_means[order(match(sub_ts_means$cluster,clusters_to_plot)),]
    
    if (num_needed_figures > 1){
      save_name<-paste0(plot_name,'_',idx,'_of_',num_needed_figures,'.svg')
    }else{
      save_name<-paste0(plot_name,'.svg')
    }
    
    sub_ts_data$labels<-factor(sub_ts_data$labels,levels = unique(sub_ts_data$labels))
    sub_ts_means$labels<-factor(sub_ts_means$labels,levels = unique(sub_ts_means$labels))
    
    cluster_num<-length(clusters_to_plot)
    number_rows<-ceiling(cluster_num/2)
    custom_height<-3*number_rows
    if (cluster_num==1){
      custom_width<-6
    }else{
      custom_width<-12
    }
    
    if(log_TP==T){
      sub_ts_data$log10_timepoint<-log10(sub_ts_data$timepoint)
      sub_ts_data$log10_timepoint[sub_ts_data$log10_timepoint=='-Inf']<-0
      sub_ts_means$log10_timepoint<-log10(sub_ts_means$timepoint)
      sub_ts_means$log10_timepoint[sub_ts_means$log10_timepoint=='-Inf']<-0
      
    }
    
    the_plot<-plot_cluster_traj_custom(object,sub_ts_data,sub_ts_means)
    svg(save_name,width=custom_width,height=custom_height)
    print(the_plot)
    dev.off()
  }
}



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


