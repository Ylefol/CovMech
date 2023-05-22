# setwd('A_Projects/EpiGen/R_Work_Folder/CovMech/')

load('results_folder/critical/TS_results/timeseries_obj_res.Rdata')

source('R/object.R')
source('R/DE_PART_results_functions.R')
source('R/MDS_GO_results_functions.R')


Elist<-readRDS('data/removed_002_127/TS_covid/processed_adjusted_cov_dta.rds')
row.names(Elist$E)=Elist$genes$GeneName
Elist<-Elist$E
healthy_patients<-colnames(Elist)[startsWith(colnames(Elist),'F')]

TS_object@sample_data<-read.csv('data/removed_002_127/TS_covid/sample_file_adjusted_w_healthy_rem_patients.csv')
TS_object@count_matrix$norm<-Elist[,TS_object@sample_data$sample]

custom_colors<-c(unname(TS_object@group_colors),'#43db0e')
names(custom_colors)<-c(names(TS_object@group_colors),'Healthy')

TS_object@group_colors<-custom_colors


create_traj_data_healthy_permutations<-function(object,custom_cmap=NULL,scale_feat=T){
  samp_dta<-object@sample_data
  h_samps<-samp_dta$sample[startsWith(samp_dta$sample,'F')]
  
  vals <- c( h_samps,h_samps,h_samps )
  vals <- unique( vals )
  perm_test<-combn( vals , 3 ) #Gives 165 possibilities, this is correct
  
  samp_df<-melt(perm_test) #Gives as a single column, gonna have to ensure sequence of three
  
  tp_vect<-c(0,48,120)
  
  healthy_samples<-data.frame(sample=samp_df$value, group=rep('Healthy',nrow(samp_df)),timepoint=rep(tp_vect,(nrow(samp_df)/length(tp_vect))))
  
  new_samp_dta<-object@sample_data[,c('sample','group','timepoint')]
  new_samp_dta<-new_samp_dta[!startsWith(new_samp_dta$sample,'F'),]
  new_samp_dta<-rbind(new_samp_dta,healthy_samples)
  
  if (is.null(custom_cmap)==T){
    my_cmap<-object@PART_results$cluster_map
  }else{
    my_cmap<-custom_cmap
  }
  
  norm_mat<-object@count_matrix$norm[row.names(my_cmap),]
  
  ts_df<-data.frame(gene_id=NULL,group=NULL,timepoint=NULL,mean_reads=NULL)
  df_list<-list()
  for (group in unique(new_samp_dta$group)){
    for (tp in unique(new_samp_dta$timepoint)){
      samples_to_merge<-new_samp_dta$sample[new_samp_dta$group==group & new_samp_dta$timepoint==tp]
      if(length(samples_to_merge)>1){
        val_vect<-rowMeans(norm_mat[,samples_to_merge])
      }else{
        val_vect<-norm_mat[,samples_to_merge,drop=F]
      }
      new_col_name<-paste0(group,'_',tp)
      df_list[[new_col_name]]<-val_vect
    }
  }
  if(scale_feat==T){
    ts_df<-sweep(data.frame(df_list),1,rowSums(data.frame(df_list)),'/')
  }else{
    ts_df<-(data.frame(df_list))
  }
  ts_df$gene_id<-row.names(ts_df)
  library(reshape)
  ts_df <- reshape2::melt(ts_df,id='gene_id')
  library(stringr)
  
  ts_df<-cbind(ts_df,str_split_fixed(ts_df$variable, "_", 2))
  # View(ts_df)
  colnames(ts_df)=c('gene_id','replicate','trans_mean','group','timepoint')
  ts_df<-ts_df[,c('gene_id','group','timepoint','trans_mean')]
  ts_df$timepoint<-as.numeric(ts_df$timepoint)
  #Add cluster association
  my_cmap$gene_id<-row.names(my_cmap)
  my_cmap<-my_cmap[,c('gene_id','cluster')]
  clust_df<-as.data.frame(table(my_cmap$cluster))
  colnames(clust_df)=c('cluster','nGenes')
  clust_df$labels<-paste0(clust_df$cluster,' – ',clust_df$nGenes,' genes |')
  
  ts_df<-merge(ts_df,my_cmap,by='gene_id')
  ts_df<-merge(ts_df,clust_df,by='cluster')
  # View(ts_df)
  ts_df$labels<-paste0(ts_df$labels,' ',ts_df$group)
  #Order the dataframe
  ts_df<-ts_df[order(ts_df$cluster,ts_df$group,ts_df$gene_id,ts_df$timepoint),]
  return(ts_df)
}

# 
# #Build custom WGCNA cmap
WGCNA_path<-'results_folder/compiled_results_16_05_2023/WGCNA_results/WGCNA_module_gene_lists/'
my_cmap<-data.frame()
for(file in list.files(WGCNA_path)){
  module_name<-strsplit(file,'.csv')[[1]]
  temp_module<-read.csv(paste0(WGCNA_path,file))
  temp_module$cluster_col<-temp_module$moduleColor
  row.names(temp_module)=temp_module$geneSymbol
  temp_module<-temp_module[,c('moduleColor','cluster_col')]
  colnames(temp_module)=c('cluster','cluster_col')

  if(nrow(my_cmap)==0){
    my_cmap<-temp_module
  }else{
    my_cmap<-rbind(my_cmap,temp_module)
  }
}
# Filter for genes which are in the normalized matrix
my_cmap<-my_cmap[row.names(my_cmap) %in% row.names(TS_object@count_matrix$norm),]
# Filter cmap for desired modules
# my_cmap<-my_cmap[my_cmap$cluster %in% c('yellow','green','black','magenta','greenyellow','midnightblue'),]

#If need to use WGCNA, input custom cmap instead of NULL
ts_data_healthy<-create_traj_data_healthy_permutations(TS_object,custom_cmap = NULL,scale_feat=T)
ts_data_healthy$group[ts_data_healthy$group=='Severe']<-'Severe'
ts_data_healthy$group[ts_data_healthy$group=='Moderate']<-'Moderate'

ts_data_healthy$labels<-gsub(x = ts_data_healthy$labels,pattern = 'Severe',replacement = 'Severe')
ts_data_healthy$labels<-gsub(x = ts_data_healthy$labels,pattern = 'Moderate',replacement = 'Moderate')

ts_data_healthy$timepoint[ts_data_healthy$timepoint==0]<-1
ts_data_healthy$timepoint[ts_data_healthy$timepoint==48]<-2
ts_data_healthy$timepoint[ts_data_healthy$timepoint==120]<-3



ts_data<-ts_data_healthy
mean_ts_data<-calculate_mean_cluster_traj(ts_data) #Calculate the mean scaled values for each cluster


clust_order<-unique(ts_data[,c('cluster','nGenes')])
clust_order<-clust_order$cluster[order(-clust_order$nGenes)]

clust_num<-length(clust_order)

num_needed_figures<-ceiling(length(clust_order)/clust_num)
#Iterate over number of necessary figures
for (idx in 1:num_needed_figures){
  #Filter the cluster_map dataframe for the required clusters
  max_clust=clust_num*idx
  if (idx==1){
    min_clust<-1
  }else{
    min_clust<-clust_num*(idx-1)+1
  }
  if (num_needed_figures > 1){
    save_name<-paste0('Ctraj_',idx,'_of_',num_needed_figures,'.svg')
  }else{
    save_name<-'Ctraj.svg'
  }
  
  
  
  clusters_to_plot<-clust_order[min_clust:max_clust]
  clusters_to_plot<-clusters_to_plot[!is.na(clusters_to_plot)]#Remove NAs
  sub_ts_data<-ts_data[ts_data$cluster %in% clusters_to_plot,]
  sub_ts_data<-sub_ts_data[order(match(sub_ts_data$cluster,clusters_to_plot)),]
  print(unique(sub_ts_data$cluster))
  
  sub_ts_means<-mean_ts_data[mean_ts_data$cluster %in% clusters_to_plot,]
  sub_ts_means<-sub_ts_means[order(match(sub_ts_means$cluster,clusters_to_plot)),]
  sub_ts_data$labels<-factor(sub_ts_data$labels,levels = unique(sub_ts_data$labels))
  sub_ts_means$labels<-factor(sub_ts_means$labels,levels = unique(sub_ts_means$labels))
  
  
  cluster_num<-length(clusters_to_plot)
  number_rows<-ceiling(cluster_num/2)
  custom_height<-3*number_rows
  if (cluster_num==1){
    custom_width<-7
  }else{
    custom_width<-14
  }
  
  mean_df<-data.frame()
  for(group in unique(sub_ts_means$group)){
    for(clust in unique(sub_ts_means$cluster)){
      sub_df<-sub_ts_means[sub_ts_means$group==group & sub_ts_means$cluster==clust,]
      new_val<-mean(sub_df$trans_mean)
      sub_df$trans_mean<-rep(new_val,nrow(sub_df))
      if(nrow(mean_df)==0){
        mean_df<-sub_df
      }else{
        mean_df<-rbind(mean_df,sub_df)
      }
    }
  }
  
  
  # plt <- ggplot(sub_ts_data, aes(y = trans_mean , x = timepoint, color = group))
  # plt <- plt + scale_color_manual(values=custom_colors) +
  #   geom_line(aes(group = gene_id), alpha = 0.4) +
  #   geom_point() +
  #   geom_line(
  #     data = sub_ts_means, lwd = 1.5, color = "grey50",
  #     aes(group = group)
  #   ) +
  #   scale_x_continuous(expand = c(0, 0)) +
  #   # scale_y_continuous(expand = c(0, 0))+
  #   facet_wrap(~labels, scales = 'free_x', ncol = 3)
  # 
  # 
  # svg('blarg.svg')
  # print(plt)
  # dev.off()

  
  
  #Remove Healthy in new sub
  no_h_dta<-sub_ts_data[sub_ts_data$group != 'Healthy',]
  no_h_mean<-sub_ts_means[sub_ts_means$group != 'Healthy',]
  
  
  # c_interest<-c('midnightblue','magenta','red','green','brown','salmon')
  c_interest<-c('C2','C3','C4','C6')
  #Order data per cluster/module
  interest<-no_h_dta[no_h_dta$cluster %in% c_interest,]
  others<-no_h_dta[!no_h_dta$cluster %in% c_interest,]
  no_h_dta<-rbind(interest,others)
  
  interest<-no_h_mean[no_h_mean$cluster %in% c_interest,]
  others<-no_h_mean[!no_h_mean$cluster %in% c_interest,]
  no_h_mean<-rbind(interest,others)
  
  healthy_means<-data.frame()
  for(label in unique(no_h_mean$labels)){
    clust=unique(no_h_mean$cluster[no_h_mean$labels==label])
    h_mean<-sub_ts_means[sub_ts_means$group == 'Healthy' & sub_ts_means$cluster==clust,]
    h_mean$trans_mean<-mean(h_mean$trans_mean)
    h_mean$labels<-rep(label,nrow(h_mean))
    if(nrow(healthy_means)==0){
      healthy_means<-h_mean
    }else{
      healthy_means<-rbind(healthy_means,h_mean)
    }
  }
  
  
  
  #Order of labels
  #Get expected order
  expect_order<-as.vector(unique(no_h_dta$labels))
  new_order<-c()
  for(i in 1:length(expect_order[c(T,F)])){
    new_order<-c(new_order,expect_order[c(F,T)][i])
    new_order<-c(new_order,expect_order[c(T,F)][i])
  }
  
  no_h_dta$cluster<-factor(no_h_dta$cluster,levels=unique(no_h_dta$cluster))
  no_h_dta$labels<-factor(no_h_dta$labels,levels=new_order)
  no_h_mean$labels<-factor(no_h_mean$labels,levels=unique(no_h_mean$labels))
  healthy_means$labels<-factor(healthy_means$labels,levels=new_order)
  

  plt <- ggplot(no_h_dta, aes(y = trans_mean , x = timepoint, color = group))
  plt <- plt +
    geom_line(aes(group = gene_id), alpha = 0.4) +
    geom_point() +
    geom_line(
      data = no_h_mean, lwd = 1.5, color = "grey50",
      aes(group = group)
    ) +
    geom_line(
      data = healthy_means, lwd = 1.5,
      aes(group = group)
    ) +
    scale_colour_manual(name = 'Groups', 
                        values =c('Severe'="#e31a1c",
                                  'Moderate'="#1f78b4",
                                  'Healthy'='#369515'), labels = c('Severe','Moderate','Healthy'))+
    # scale_x_continuous(expand = c(0, 0)) + 
    scale_x_continuous(name ="timepoints",
                    breaks=c(1,2,3),
                    labels=c("T1", "T2", "T3"))+
    # scale_x_discrete(breaks=c('0','48','168'),
    #                    labels=c("Day 1", "Day 3", "Day 8"))+
    ylab("Scaled expression") +
    facet_wrap(~labels, scales = 'free_x', ncol = 4)+
    theme(strip.text.x = element_text(size = 12,face="bold"))+
    theme(axis.text = element_text(size = 16))+
    theme(legend.position = "none") +
    theme(axis.title=element_text(size=14,face="bold")) +
    theme(axis.title.x = element_blank())
  
  svg(save_name,width = custom_width,height=custom_height)
  print(plt)
  dev.off()

}



