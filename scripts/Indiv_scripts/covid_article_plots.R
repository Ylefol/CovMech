#Custom figures for CoV-Mech article
# setwd('~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/')
#Dotplots for figure 1D/2 – merger of timepoints
# load('results_folder/critical/cluster_profiler_results/severe.cov.T1_vs_moderate.cov.T1/cluster_profiler_objects.RData')
load('results_folder/critical/cluster_profiler_results/all.cov.T1_vs_healthy.controls///cluster_profiler_objects.RData')
# load('results_folder/viremia/cluster_profiler_results/viremia.cov.D1_vs_non_viremia.cov.D1/cluster_profiler_objects.RData')

TP1_gse<-gse
# load('results_folder/critical/cluster_profiler_results/severe.cov.T2_vs_moderate.cov.T2/cluster_profiler_objects.RData')
load('results_folder/critical/cluster_profiler_results/all.cov.T2_vs_healthy.controls////cluster_profiler_objects.RData')
# load('results_folder/viremia/cluster_profiler_results/viremia.cov.D3_vs_non_viremia.cov.D3/cluster_profiler_objects.RData')

TP2_gse<-gse
# load('results_folder/critical/cluster_profiler_results/severe.cov.T3_vs_moderate.cov.T3//cluster_profiler_objects.RData')
load('results_folder/critical/cluster_profiler_results/all.cov.T3_vs_healthy.controls////cluster_profiler_objects.RData')
# load('results_folder/viremia/cluster_profiler_results/viremia.cov.D8_vs_non_viremia.cov.D8/cluster_profiler_objects.RData')

TP3_gse<-gse

rm(gene_list,go_enrich_down,go_enrich_up,gse)
source('R/my_ClusterProfiler_functions.R')

custom_GO_selection<-c('neutrophil mediated immunity','neutrophil degranulation','type I interferon signaling pathway',
                       'response to type I interferon','response to interferon-gamma','interferon-gamma-mediated signaling pathway',
                       'defense response to virus','cytokine-mediated signaling pathway','cytokine production',
                       'wound healing','coagulation','platelet degranulation','platelet activation','humoral immune response',
                       'RNA catabolic process','translation','ribosome biogenesis','antigen processing and presentation of peptide antigen via MHC class II',
                       'MHC class II protein complex assembly','T cell receptor signaling pathway')
#Need to build a df with everything, use rbind to slap them together
#at their different timepoints.
#Start by getting the top 10 of both positively and negatively enriched

# gse_neg<-TP1_gse@result$Description[TP1_gse@result$enrichmentScore<0][1:20]
# gse_pos<-TP1_gse@result$Description[TP1_gse@result$enrichmentScore>0][1:20]
# TP1_desc<-c(gse_neg[1:10],gse_pos[1:10])
TP1_res<-TP1_gse@result[TP1_gse@result$Description %in% custom_GO_selection,c('ID','Description','setSize','enrichmentScore',"p.adjust")]
TP1_res[['-log10(padj)']]<--log10(TP1_res$p.adjust)
row.names(TP1_res)=c(1:nrow(TP1_res))
TP1_res$timepoints<-rep(x = 'TP1',nrow(TP1_res))

# gse_neg<-TP2_gse@result$Description[TP2_gse@result$enrichmentScore<0][1:20]
# gse_pos<-TP2_gse@result$Description[TP2_gse@result$enrichmentScore>0][1:20]
# TP2_desc<-c(gse_neg[1:10],gse_pos[1:10])
TP2_res<-TP2_gse@result[TP2_gse@result$Description %in% custom_GO_selection,c('ID','Description','setSize','enrichmentScore',"p.adjust")]
TP2_res[['-log10(padj)']]<--log10(TP2_res$p.adjust)
row.names(TP2_res)=c(nrow(TP1_res)+1:nrow(TP2_res))
TP2_res$timepoints<-rep(x = 'TP2',nrow(TP2_res))


# gse_neg<-TP3_gse@result$Description[TP3_gse@result$enrichmentScore<0][1:20]
# gse_pos<-TP3_gse@result$Description[TP3_gse@result$enrichmentScore>0][1:20]
# TP3_desc<-c(gse_neg[1:10],gse_pos[1:10])
TP3_res<-TP3_gse@result[TP3_gse@result$Description %in% custom_GO_selection,c('ID','Description','setSize','enrichmentScore',"p.adjust")]
TP3_res[['-log10(padj)']]<--log10(TP3_res$p.adjust)
row.names(TP3_res)=c(nrow(TP2_res)+nrow(TP1_res)+1:nrow(TP3_res))
TP3_res$timepoints<-rep(x = 'TP3',nrow(TP3_res))

full_data<-rbind(TP1_res,TP2_res,TP3_res)

#Set-up necessary functions
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}


label_format<-20
label_func <- default_labeller(label_format)
if(is.function(label_format)) {
  label_func <- label_format
}

# full_data$timepoint[full_data$timepoint=='TP1']<-'D1'
# full_data$timepoint[full_data$timepoint=='TP2']<-'D3'
# full_data$timepoint[full_data$timepoint=='TP3']<-'D8'

full_data$Description<-factor(full_data$Description,levels=rev(custom_GO_selection))

#With the full data, gotta make the plot now.
plt<-ggplot(full_data, aes(x=timepoints, y=Description, size=`-log10(padj)`, color=enrichmentScore)) +
  geom_point() +
  scale_fill_gradientn(colours=rev(c("red","white","blue")),
                       # values=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                       breaks=c(-1,0,1),
                       # breaks=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                       limits=c(-1, 1), oob=squish,
                       guide=guide_colorbar(reverse=F), name = 'enrichmentScore',aesthetics = "colour")+ 
  # scale_color_continuous(low="blue", high="gray", name = color,
  #                        limits=c(-1, 0),
  #                        guide='legend') +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 50)) +
  ylab(NULL) + ggtitle('Critical vs Mopderate') + theme_dose(15) +
  scale_size(range=c(3, 8))

ggsave(file='merged_tp_dotplot.png',plt,width = 10, height = 8,dpi=300)




#Plot original C6 genes for the bad mech vent group

# setwd("~/A_Projects/EpiGen/R_Work_Folder/TS_dev_zone")
# 
# og_cmap<-read.csv('Results_folder/covid_results/TS_cov_results/PART_results/PART_heat_cmap.csv',row.names=1)
# new_cmap<-og_cmap[og_cmap$cluster=='C6',]
# 
# load('Results_folder/covid_results/TS_results_rem_bad_mechvent/timeseries_obj_res.Rdata')
# source('R/DE_PART_results_functions.R')
# my_ts_data<-calculate_cluster_traj_data(TS_object,custom_cmap = new_cmap,scale_feat=T)
# my_mean_ts_data<-calculate_mean_cluster_traj(my_ts_data)
# the_plot<-plot_cluster_traj(TS_object,my_ts_data,my_mean_ts_data)
# ggsave(file='C6_bad_mech_vent.png',the_plot)






#To get the highres PART heatmap, I should be able to simply source the code, load the object and run it?
#Adjusting the dpi may a bit more complicated




#custom PART plot function
source('R/DE_PART_results_functions.R')
PART_heat_map_custom<-function(object, heat_name='custom_heat_map'){
  
  clust_ordered <- unique(as.character(object@PART_results$cluster_info[['calculated_clusters']]$lab.hatK[object@PART_results$cluster_info[['clustered_rows']]$order]))
  #Cluster illustration stored as rowAnnotation
  row_annot <- ComplexHeatmap::rowAnnotation(gene_cluster = as.character(object@PART_results$cluster_info[['calculated_clusters']]$lab.hatK), 
                                             col = list(gene_cluster=object@PART_results$cluster_info[['colored_clust_rows']]),
                                             show_annotation_name=F,
                                             annotation_legend_param = list(title = "clusters", at = clust_ordered,
                                                                            labels = unique(object@PART_results$cluster_map$cluster)))
  
  #Create top annotations
  top_annot_results<-prepare_top_annotation_PART_heat(object)
  top_annot_labels<-top_annot_results[[1]]
  top_annot_no_labels<-top_annot_results[[2]]
  col_split<-top_annot_results[[3]]
  
  #Prepare a 'gap_vect' to split the different elements of the heatmap
  col_split_vect<-levels(col_split)
  gap_vect<-c()
  for (idx in 1:length(col_split_vect)){
    if (idx!=length(col_split_vect)){
      val1<-unlist(strsplit(as.character(col_split_vect[idx]),'_'))[1]
      val2<-unlist(strsplit(as.character(col_split_vect[idx+1]),'_'))[1]
      if (val1!=val2){
        gap_vect<-c(gap_vect,4)
      }else{
        gap_vect<-c(gap_vect,0.5)
      }
    }
  }
  
  #Create the legend for the groups
  lgd = Legend(labels = names(object@group_colors), title = "groups", 
               legend_gp = gpar(fill = unname(object@group_colors)))

  #Plot the heatmap
  reso <- 300
  length <- 3.25*reso/72
  width<-length-(length/4)
  
  png(paste0(heat_name,".png"),units='in',width=width,height=length,res = reso)
  draw(
    ComplexHeatmap::Heatmap(
      object@PART_results$part_matrix, name = "Z-score", cluster_columns = F,
      cluster_rows=object@PART_results$cluster_info[['clustered_rows']], show_column_dend = T,
      row_names_gp = grid::gpar(fontsize = 8),left_annotation = row_annot,
      show_row_names = F,top_annotation = top_annot_no_labels,column_split = col_split,cluster_column_slices = T,
      column_gap = unit(gap_vect, "mm"),show_column_names = F,border=F,column_title = NULL),
    annotation_legend_list = lgd
  )
  trash<-capture.output(dev.off())#Capture output to prevent print
}

load('TS_results/timeseries_obj_res.Rdata')


PART_heat_map_custom(TS_object)








# #Get the GOs for TS, append into excel file as sheets
# #Only keep relevant information
# library("xlsx")
# file_locations<-'TS_results/gprofiler_results/data_files/'
# excel_file_name<-'TS_GOs.xlsx'
# 
# for(item in list.files(file_locations)){
#   name<-strsplit(item,'_data.csv')[[1]]
#   print(name)
#   my_dta<-read.csv(paste0(file_locations,item),row.names = 1)
#   my_dta<-my_dta[,c('term_id','term_name','p_value','term_size','intersection_size')]
#   
#   if(excel_file_name %in% list.files()==F){
#     write.xlsx(my_dta, excel_file_name, sheetName = name, 
#                col.names = T, row.names = T, append = F)
#   }else{
#     write.xlsx(my_dta, excel_file_name, sheetName = name, 
#                col.names = T, row.names = T, append = T)
#   }
# }

# 
# 
# library("xlsx")
# file_locations<-'WGCNA_results/gprofiler_results/data_files/'
# excel_file_name<-'WGCNA_GOs.xlsx'
# 
# for(item in list.files(file_locations)){
#   name<-strsplit(item,'_data.csv')[[1]]
#   print(name)
#   my_raw_dta<-read.csv(paste0(file_locations,item))
#   my_dta<-my_raw_dta[my_raw_dta$source=='GO:BP',]
#   my_dta<-my_dta[,c('term_id','term_name','p_value','term_size','intersection_size')]
#   
#   if(nrow(my_dta)>0){
#     if(excel_file_name %in% list.files()==F){
#       write.xlsx(my_dta, excel_file_name, sheetName = name, 
#                  col.names = T, row.names = F, append = F)
#     }else{
#       write.xlsx(my_dta, excel_file_name, sheetName = name, 
#                  col.names = T, row.names = , append = T)
#     }
#   }
# 
# }