setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")
# source('code/functions/my_ClusterProfiler_functions.R')
source('R/DE_PART_results_functions.R')

my_list<-c('all.cov.D1_vs_healthy.controls','all.cov.D3_vs_healthy.controls','all.cov.D8_vs_healthy.controls')
sig_file<-'DE_significant_1_005.csv'
path<-'DE_results/'

prep_object<-function(my_list,sig_file,path){
  DEG_list<-list()
  for_replicate_vect<-c()
  for(exp in my_list){
    timepoint<-strsplit(strsplit(exp,'_vs_')[[1]][1],'\\.')[[1]][3]
    DEG_list[[exp]]<-read.csv(paste0(path,exp,'/',sig_file),check.names = F)
    found_replicates<-colnames(DEG_list[[exp]])[7:length(colnames(DEG_list[[exp]]))]
    print(found_replicates)
    new_replicates<-c()
    for (rep in found_replicates){
      rep<-strsplit(rep,'_')[[1]][1]
      rep<-paste0(rep,'_',timepoint)
      new_replicates<-c(new_replicates,rep)
    }
    
    colnames(DEG_list[[exp]])[7:length(colnames(DEG_list[[exp]]))]<-new_replicates
    # found_replicates<-colnames(DEG_list[[exp]])[8:length(colnames(DEG_list[[exp]]))]
    for_replicate_vect<-c(for_replicate_vect,new_replicates)
    
  }
  rep_vect<-c()
  group_vect<-c()
  for (samp in for_replicate_vect){
    rep_vect<-c(rep_vect,strsplit(samp,'_')[[1]][1])
    if(startsWith(samp,'F')==T){
      group_vect<-c(group_vect,'Healthy')
    }else{
      group_vect<-c(group_vect,'Covid')
    }
  }
  group_df<-data.frame(Group=group_vect,Replicate=rep_vect)
  # row.names(group_df)=group_df$Replicate
  
  names(rep_vect)<-for_replicate_vect
  
  return(list(DEG_list,group_df,rep_vect))
}


create_heatmap_matrix_temp<-function(rep_vect,group_names,group_df,DEG_list){
  #Create empty list
  heat_list<-list()
  
  #Establish name of experiments to search for
  my_groups<-paste0(group_names[2],'_vs_',group_names[1])
  
  
  gene_vector<-c()
  #Create empty df and vector
  heat_df<-as.data.frame(NULL)
  l2fc_vect<-c()
  
  
  exps_interest<-names(DEG_list)
  not_interest<-c('merged_df','sorted_df','unique_genes')
  exps_interest<-exps_interest[!exps_interest %in% not_interest]
  
  list_of_formatted_exp<-list()
  #Iterate over experiments/timepoints
  for (tp_group in exps_interest){
    # tp_group<-paste0(my_groups,'_TP_',tp)#Create name to look for
    temp_df<-as.data.frame(DEG_list[[tp_group]])#Extract df associated to timepoint
    
    genes_mod<-temp_df$gene_id#Extract gene names
    genes_mod<-paste0(genes_mod,'_',tp_group)#Modify gene names by adding time point
    
    #Extract and name log2FoldChange
    temp_l2fc<-temp_df$log2FoldChange
    names(temp_l2fc)<-genes_mod
    
    temp_df<-temp_df[,7:ncol(temp_df)]#Select for experiment related data
    
    row.names(temp_df)=genes_mod
    temp_df<-t(temp_df)#transpose
    row.names(temp_df)=unname(rep_vect[row.names(temp_df)])#rename using replicates
    
    
    list_of_formatted_exp[[tp_group]]<-temp_df
    
    l2fc_vect<-c(l2fc_vect,temp_l2fc)
    temp_gene_vector<-rep(paste0(tp_group,' (',length(genes_mod),')'),length(genes_mod))
    gene_vector<-c(gene_vector,temp_gene_vector)
  }
  
  all_rows<-c()
  for(exp in names(list_of_formatted_exp)){
    all_rows<-c(all_rows,row.names(list_of_formatted_exp[[exp]]))
  }
  all_rows<-unique(all_rows)
  heat_df<-data.frame(NULL)
  #Start by building each dataframe individually, then c bind them since the rows will be the same
  for (exp in names(list_of_formatted_exp)){
    missing_rows<-all_rows[!all_rows %in% row.names(list_of_formatted_exp[[exp]])]
    if (length(missing_rows)>0){
      NA_vect<-rep(NA,length(colnames(list_of_formatted_exp[[exp]])))
      NA_df<-data.frame(NA_vect)
      NA_df<-do.call("cbind", replicate(length(missing_rows), NA_df, simplify = FALSE))
      NA_df<-t(NA_df)
      row.names(NA_df)<-missing_rows
      list_of_formatted_exp[[exp]]<-rbind(list_of_formatted_exp[[exp]],NA_df)
    }
    list_of_formatted_exp[[exp]]<-list_of_formatted_exp[[exp]][all_rows,]
    
    if (nrow(heat_df)==0){
      heat_df<-list_of_formatted_exp[[exp]]
    }else{
      heat_df<-cbind(heat_df,list_of_formatted_exp[[exp]])
    }
  }
  # my_group_df<-as.data.frame(TimeSeriesExperiment::groups(time_object))
  # my_group_df<-cbind(my_group_df,unname(replicates(time_object)))
  # colnames(my_group_df)<-c('Group','Replicate')
  group_df<-unique(group_df)
  row.names(group_df)=group_df$Replicate
  
  group_vector<-group_df[row.names(heat_df),]
  print(group_vector)
  group_vector<-group_vector$Group
  print(group_vector)
  #Fill and return result list
  heat_list[['main_matrix']]<-heat_df
  heat_list[['l2fc_vector']]<-l2fc_vect
  heat_list[['group_vector']]<-group_vector
  heat_list[['gene_vector']]<-gene_vector
  return(heat_list)
}

plot_custom_DE_covid_heat <-function(heat_mat,col_split,row_splits,l2fc_col, log_transform,
                                  save_path='',plot_file_name='custom_heatmap',
                                  custom_width=15,custom_height=5){
  

  bottom_histo = HeatmapAnnotation('log10(FC)' = anno_barplot(l2fc_col))
  count_legend<-'log10(intensity)'
  
  
  
  fill_set_regions=gpar(fill=2:(length(unique(col_split))+1))
  top_annot = HeatmapAnnotation(foo = anno_block(gp = fill_set_regions, labels = unique(col_split)))
  
  fill_set_groups=gpar(fill = c("#fdc086", "#beaed4"))

  left_annot = rowAnnotation(foo = anno_block(gp = fill_set_groups, labels = unique(row_splits)))

  
  legend_label<-unique(col_split)
  #Calculate ideal distance between legend points in mm, the max distance will be used
  gap_vect<-nchar(levels(legend_label))*1.8
  
  #Create legend
  lgd = Legend(labels = legend_label, title = "Experiments", legend_gp = fill_set_regions,nrow = 1, 
               title_position = "leftcenter",gap = unit(max(gap_vect), "mm"))
  
  
  save_name_svg<-paste0(save_path,plot_file_name,'.svg')
  
  svg(save_name_svg,width=custom_width,height=custom_height)
  draw(Heatmap(heat_mat,name=count_legend,cluster_rows = T,cluster_columns = T,
               show_column_names = F,show_row_names = F,row_names_side='left',
               row_split = row_splits,column_split=col_split,
               border=T,na_col = 'gray',column_title = NULL,row_title = NULL,
               top_annotation = top_annot,left_annotation=left_annot,
               bottom_annotation = bottom_histo,
               cluster_column_slices=F,cluster_row_slices = F)#,
       # annotation_legend_list = lgd,
       # annotation_legend_side = 'bottom'
  )
  dev.off()
  
  save_name_png<-paste0(save_path,plot_file_name,'.png')
  
  #width and height *96 to convert inches to pixels
  png(save_name_png,width=custom_width*96,height=custom_height*96)
  draw(Heatmap(heat_mat,name=count_legend,cluster_rows = T,cluster_columns = T,
               show_column_names = F,show_row_names = F,row_names_side='left',
               row_split = row_splits,column_split=col_split,
               border=T,na_col = 'gray',column_title = NULL,row_title = NULL,
               top_annotation = top_annot,left_annotation=left_annot,
               bottom_annotation = bottom_histo,
               cluster_column_slices=F,cluster_row_slices = F)#,
       # annotation_legend_list = lgd,
       # annotation_legend_side = 'bottom'
  )
  dev.off()
  
  
  return(row_order(draw(Heatmap(heat_mat,name=count_legend,cluster_rows = T,cluster_columns = T,
                                show_column_names = F,show_row_names = F,row_names_side='left',
                                row_split = row_splits,column_split=col_split,
                                border=T,na_col = 'gray',column_title = NULL,row_title = NULL,
                                top_annotation = top_annot,left_annotation=left_annot,
                                bottom_annotation = bottom_histo,
                                cluster_column_slices=F,cluster_row_slices = F)#,
                        # annotation_legend_list = lgd,
                        # annotation_legend_side = 'bottom'
  )))
  
}



my_objs<-prep_object(my_list,sig_file,path)

bloop<-create_heatmap_matrix_temp(my_objs[[3]],c('Healthy','Covid'),my_objs[[2]],my_objs[[1]])
mat_cols<-colnames(bloop[[1]])
mat_cols<-gsub(x = mat_cols,pattern = '_all.cov.D1_vs_healthy.controls',replacement = ' Day 1')
mat_cols<-gsub(x = mat_cols,pattern = '_all.cov.D3_vs_healthy.controls',replacement = ' Day 3')
mat_cols<-gsub(x = mat_cols,pattern = '_all.cov.D8_vs_healthy.controls',replacement = ' Day 8')
colnames(bloop[[1]])=mat_cols


#Now take care of 2 in the same way
to_sub<-names(bloop[[2]])
to_sub<-gsub(x = to_sub,pattern = '_all.cov.D1_vs_healthy.controls',replacement = ' Day 1')
to_sub<-gsub(x = to_sub,pattern = '_all.cov.D3_vs_healthy.controls',replacement = ' Day 3')
to_sub<-gsub(x = to_sub,pattern = '_all.cov.D8_vs_healthy.controls',replacement = ' Day 8')
names(bloop[[2]])=to_sub


#Now number 4
to_sub<-bloop[[4]]
to_sub[to_sub=='all.cov.D1_vs_healthy.controls (1948)']<-'Day 1 (1948)'
to_sub[to_sub=='all.cov.D3_vs_healthy.controls (1621)']<-'Day 3 (1621)'
to_sub[to_sub=='all.cov.D8_vs_healthy.controls (1044)']<-'Day 8 (1044)'
bloop[[4]]<-to_sub

log_transform=T
plot_file_name='covid_vs_healthy'
temp_list<-prepare_heat_data(bloop,log_transform = log_transform)

my_heat_mat<-temp_list[['heat_matrix']]
my_region_split<-temp_list[['region_split']]
my_group_split<-temp_list[['group_split']]
my_l2fc_vect<-temp_list[['l2fc_vector']]

my_l2fc_vect<-log_transform_l2fc_vect(my_l2fc_vect)

sample_order<-plot_custom_DE_covid_heat(my_heat_mat,my_region_split,my_group_split,my_l2fc_vect,log_transform = log_transform,
                                        save_path='',plot_file_name = plot_file_name)



row_order<-row.names(my_heat_mat)[c(sample_order$Covid,sample_order$Healthy)]

samp_dta<-read.csv('data/TS_covid/sample_file_no_healthy_rem_patients.csv')
samp_dta$sample<-unlist(strsplit(samp_dta$sample,'_'))[c(T,F)]
samp_dta<-unique(samp_dta[,c('sample','group')])

health_samps<-row.names(my_heat_mat)[!row.names(my_heat_mat) %in% samp_dta$sample]
healthy_df<-data.frame(sample=health_samps,group=rep('Healthy',length(health_samps)))

full_samp_dta<-rbind(samp_dta,healthy_df)
row.names(full_samp_dta)=full_samp_dta$sample

full_samp_dta<-full_samp_dta[row_order,]

write.csv(file = 'heat_order.csv',full_samp_dta)
