#Rmarkdown functions

find_most_variable_cluster<-function(object,mean_ts_data){
  
  #Retrieve clusters which have gprofiler results
  clusts_interest<-names(object@Gprofiler_results)
  mean_ts_data<-mean_ts_data[mean_ts_data$cluster %in% clusts_interest,]
  
  dif_clust<-c()
  for(clust in unique(mean_ts_data$cluster)){
    dif_vect<-c()
    for(tp in unique(mean_ts_data$timepoint)){
      val_exp<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==TS_object@group_names[1]]
      val_control<-mean_ts_data$trans_mean[mean_ts_data$timepoint==tp & mean_ts_data$cluster==clust & mean_ts_data$group==TS_object@group_names[2]]
      
      dif_val<-abs(val_exp-val_control)
      dif_vect<-c(dif_vect,dif_val)
    }
    dif_clust_val<-sum(dif_vect)
    dif_clust<-c(dif_clust,dif_clust_val)
  }
  names(dif_clust)=unique(mean_ts_data$cluster)
  target_clust<-names(dif_clust)[unname(dif_clust)==max(dif_clust)]
  
  return(target_clust)
  
}


create_DEG_df<-function(time_object){
  #Get DEGs
  DEG_amount<-c()
  names_exp<-c()
  gene_vect<-c()
  for(DE_type in names(time_object@DE_results)){
    for(exp in names(time_object@DE_results[[DE_type]])){
      gene_vect<-c(gene_vect,time_object@DE_results[[DE_type]][[exp]][['DE_sig_data']]$gene_id)
      DEG_amount<-c(DEG_amount,nrow(time_object@DE_results[[DE_type]][[exp]][['DE_sig_data']]))
      exp<-paste0(exp,' (',DE_type,')')
      names_exp<-c(names_exp,exp)
    }
  }
  DEG_amount<-c(DEG_amount,paste0('**',length(unique(gene_vect)),'**'))
  names_exp<-c(names_exp,'**total unique genes**')
  DEG_summary<-data.frame(`experiment name`=names_exp,`number of DEGs`=DEG_amount)
  
  return(DEG_summary)
}
                        