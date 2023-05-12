
library(DESeq2)
library(ggplot2)
library(ggrepel)

# wrapper functions -------------------

#' Wrapper function which plots the trajectory of WGCNA modules for a single time
#' series object. The functions retrieves the WGCNA module genes, and calls the necessary
#' functions to process, format, and plot the trajectory data.
#'
#' @param path_to_ts string to the path of the time series object
#' @param path_to_modules string for the path to the folder containing the module
#' files of the WGCNA analysis
#' @param file_name The name to be given to the file(s) that will be produced
#' by the function
#'
#' @export
wrapper_WGCNA_traj<-function(path_to_ts, path_to_modules, file_name){
  load(file = path_to_ts)
  module_cmap<-data.frame(NULL)
  for (module in list.files(path_mods)){
    mod_res<-read.csv(paste0(path_mods,'/',module),row.names = 1)
    row.names(mod_res)=mod_res$geneSymbol
    mod_res$cluster<-rep(strsplit(module,'.csv')[[1]][1],nrow(mod_res))

    mod_res<-mod_res[,c('cluster','moduleColor')]
    colnames(mod_res)=c('cluster','cluster_col')

    if(nrow(module_cmap)==0){
      module_cmap<-mod_res
    }else{
      module_cmap<-rbind(module_cmap,mod_res)
    }
  }

  module_cmap<-module_cmap[row.names(module_cmap) %in% row.names(TS_object@count_matrix$norm),]


  my_ts_data<-calculate_cluster_traj_data(TS_object,custom_cmap=module_cmap,scale_feat=T)
  my_mean_ts_data<-calculate_mean_cluster_traj(my_ts_data)

  wrapper_cluster_trajectory(TS_object,my_ts_data,my_mean_ts_data,plot_name = file_name)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Function which creates a dotplot of found children to specified GO ancestors
#' for merged timeseries objects
#'
#' The function loads in the necessary time series objects and searches for any
#' children of the specified ancestors. It then plots these children based on the
#' clusters they were found in.
#'
#' @param vect_of_ts_paths a vector of paths to different time series objects
#' @param target_ancestors A vector of GO IDs representing the ancestors being queried
#' @param ontology The ontology to search in gprofiler format (ex: GO:BP)
#' @param plot_title The title that will be given to the dotplot
#'
#' @export
merged_ancestor_dotplot<-function(vect_of_ts_paths,target_ancestors,ontology='GO:BP',plot_title=''){
  my_ancestor_df<-data.frame(NULL)
  for(TS in vect_of_ts_paths){
    load(TS)
    name_group<-TS_object@group_names[1]
    GO_clusters<-gprofiler_cluster_analysis(TS_object,ontology,save_path = NULL)[['GO_df']]
    onto_ancestor<-unlist(str_split(ontology,':'))[2]
    temp_ancestor<-find_relation_to_ancestors(target_ancestors,GO_clusters,ontology = onto_ancestor)

    if(nrow(my_ancestor_df)==0){
      my_ancestor_df<-temp_ancestor
    }else{
      my_ancestor_df<-rbind(my_ancestor_df,temp_ancestor)
      my_ancestor_df$group_name<-factor(my_ancestor_df$group_name,levels=unique(my_ancestor_df$group_name))
    }
  }

  plt <-ggplot(my_ancestor_df, aes(x = group_name, y = term_name, color=`-log10(padj)`)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_point(aes(size=term_size),show.legend = T) +
    scale_size(name   = "GO size",
               breaks = c(min(my_ancestor_df$term_size),max(my_ancestor_df$term_size)),
               labels = c(min(my_ancestor_df$term_size),max(my_ancestor_df$term_size))) +
    # scale_size(name   = "GO size",
    #            breaks = c(min(ancestor_df$term_size),max(ancestor_df$term_size)),
    #            labels = c(min(ancestor_df$term_size),max(ancestor_df$term_size))) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , "\n"),width = 60))+
    theme_bw(base_size=16) +
    theme(axis.text.y=element_text(size = 20)) +
    theme(axis.text.x=element_text(angle=-45,vjust=0,size=20)) +
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=25), #change legend title font size
          legend.text = element_text(size=25)) +#change legend text font size
    ylab("") +
    xlab("") +
    ggtitle(plot_title)#+
  # theme(legend.position = "none")

  return_list<-list(ancestor_plot=plt,ancestor_df=my_ancestor_df)
  return(return_list)
}

# processing/formatting functions -------------------


#' Function which loads several time series object and merges them to create
#' multi group objects. For downstream analyses, the function returns raw and normalized
#' count matrices, the deseq2 object, and a dataframe of sample data for all samples
#' contained in the objects
#'
#' The function only retrieves elements from RNAseq based itmeseries analyses as
#' it utilizes DESeq2
#'
#' @param vect_of_ts_paths a vector of paths to different time series objects
#' @param group_names a vector for the names of the groups contained in the objects
#' that should be merged together
#'
#' @return a list of merged objects
#'
#' @export
load_several_RNAseq_groups<-function(vect_of_ts_paths,group_names){
  raw_matrix<-data.frame(NULL)
  samp_dta<-data.frame(NULL)
  for(TS in vect_of_ts_paths){
    load(TS)
    if(nrow(raw_matrix)==0){
      raw_matrix<-TS_object@count_matrix$raw
      samp_dta<-TS_object@sample_data
    }else{
      missing<-TS_object@sample_data[!TS_object@sample_data$sample %in% samp_dta$sample,]
      samp_dta<-rbind(samp_dta,missing)
      raw_matrix<-cbind(raw_matrix,TS_object@count_matrix$raw[,missing$sample])
    }
  }
  #Assign Condition to all of the columns/samples/patients
  my_condition<-c()
  for(group in group_names){
    my_condition<-c(my_condition,rep(group,nrow(samp_dta[samp_dta$group==group,])))
  }
  condition<-factor(my_condition)


  #Create a coldata frame and instantiate the DESeqDataSet
  col_data <- data.frame(row.names=colnames(raw_matrix),condition)
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(raw_matrix), colData=col_data, design=~condition)
  dds = estimateSizeFactors(dds)

  norm_counts <- as.data.frame(counts(dds,normalized=TRUE))
  norm_counts<- na.omit(norm_counts)

  return_list<-list(dds_obj=dds,raw_counts=raw_matrix,normalized_counts=norm_counts,sample_data=samp_dta)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which calculates the trajectory data for a single gene using different
#' time series objects and groups
#'
#' @param vect_of_ts_paths a vector of paths to different time series objects
#' @param gene_name The name of a gene that is in the normalized matrix of all
#' time series objects submitted
#'
#' @return The trajectorÿ data in plotting format for the submitted gene
#'
#' @export
data_merged_gene_traj<-function(vect_of_ts_paths,gene_name){
  mean_data<-data.frame(NULL)
  for (TS in vect_of_ts_paths){
    load(TS)
    #Create a list of genes of interest
    my_dta<-TS_object@count_matrix$norm[gene_name,]
    library(reshape2)
    my_dta<-melt(my_dta)
    colnames(my_dta)=c('sample','reads')
    my_dta<-merge(my_dta,TS_object@sample_data,by='sample')

    exp<-TS_object@group_names[1]
    mean_data_list<-list()
    for (group in unique(my_dta$group)){
      for (tp in unique(my_dta$timepoint)){
        new_val=mean(my_dta$reads[my_dta$group==group & my_dta$timepoint==tp])
        mean_data_list[['reads']]<-c(mean_data_list[['reads']],new_val)
        mean_data_list[['timepoint']]<-c(mean_data_list[['timepoint']],tp)
        mean_data_list[['label']]<-c(mean_data_list[['label']],exp)
        if (group==exp){
          mean_data_list[['group']]<-c(mean_data_list[['group']],'Activated')
        }else{
          mean_data_list[['group']]<-c(mean_data_list[['group']],'Control')
        }
      }
    }
    if (nrow(mean_data)==0){
      mean_data<-as.data.frame(mean_data_list)
    }else{
      mean_data<-rbind(mean_data,as.data.frame(mean_data_list))
    }
  }
  return(mean_data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates a dataframe containing the gene names of the clusters of
#' interest. The functions takes in a vector of time series objects and a vector of
#' clusters of interest. Both vectors must be the same size and the clusters of interest
#' must be in the order of the given time series paths. The genes for the first cluster
#' will be retrieved from the first path, the second cluster from the second
#' path etc... This means that there can only be one cluster of interest per timeseries
#' object.
#' All genes will be stored under the umbrella of 'C1'
#'
#' @param vect_of_ts_paths a vector of paths to different time series objects
#' @param clusters_interest a vector of cluster names
#'
#' @return a dataframe with the gene names and cluster name ('C1')
custom_cmap_diff_groups<-function(vect_of_ts_paths,clusters_interest){
  genes_interest<-c()
  idx_iter=1
  for(TS in vect_of_ts_paths){
    load(TS)
    genes_interest<-c(genes_interest,
                      row.names(TS_object@PART_results$cluster_map[TS_object@PART_results$cluster_map$cluster==clusters_interest[idx_iter],]))
    idx_iter=idx_iter+1
  }
  genes_interest<-unique(genes_interest)

  cmap<-data.frame(gene_id=genes_interest,cluster=rep('C1',length(genes_interest)))
  row.names(cmap)=cmap$gene_id
  return(cmap)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which retrieves all the genes of clusters belonging to different groups
#' and prepares a dataframe to be able to plot those genes in the different groups given
#' The preparation involves some dataframe formating as well as the scaling of the
#' data
#'
#' @param RNAseq_dta_objects List of merged objects as created by \code{load_several_RNAseq_groups}
#' @param cmap a dataframe indicating which clusters are of interest. Often produced
#' by \code {custom_cmap_diff_groups}
#'
#' @return The scaled and formatted dataframe
#'
#' @export
cluster_traj_merged_groups<-function(RNAseq_dta_objects,cmap){
  #Variable setup
  my_samp_dta<-RNAseq_dta_objects[['sample_data']]
  norm_mat<-RNAseq_dta_objects[['normalized_counts']]
  my_cmap=cmap

  ts_df<-data.frame(gene_id=NULL,group=NULL,timepoint=NULL,mean_reads=NULL)
  df_list<-list()
  for (group in unique(my_samp_dta$group)){
    for (tp in unique(my_samp_dta$timepoint)){
      samples_to_merge<-my_samp_dta$sample[my_samp_dta$group==group & my_samp_dta$timepoint==tp]
      val_vect<-rowMeans(norm_mat[,samples_to_merge])
      new_col_name<-paste0(group,'_',tp)
      df_list[[new_col_name]]<-val_vect
    }
  }
  ts_df<-sweep(data.frame(df_list),1,rowSums(data.frame(df_list)),'/')


  ts_df$gene_id<-row.names(ts_df)
  library(reshape)
  ts_df <- reshape2::melt(ts_df,id='gene_id')
  library(stringr)
  ts_df<-cbind(ts_df,str_split_fixed(ts_df$variable, "_", 2))
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
  ts_df$labels<-paste0(ts_df$labels,' ',ts_df$group)
  #Order the dataframe
  ts_df<-ts_df[order(ts_df$cluster,ts_df$group,ts_df$gene_id,ts_df$timepoint),]
  ts_df$gene_id<-paste0(ts_df$gene_id,ts_df$group)

  return(ts_df)
}



# plotting functions -------------------

#' Function which takes in merged time series objects from the \code{load_several_RNAseq_groups}
#' function, and then plots a pca plot for all those samples.
#'
#' The PCA plot is done using the vst funciton of DESeq2
#'
#' @param RNAseq_obj Object created by \code{load_several_RNAseq_groups}
#' @param group_order Vector of group names in the order in which they should appear
#' in the legend (from top to bottom)
#'
#' @return The ggplot2 object for the pca plot
#'
#' @export
plot_group_PCA_RNAseq<-function(RNAseq_obj,group_order){

  dds<-RNAseq_obj[['dds_obj']]
  my_samp_dta<-RNAseq_obj[['sample_data']]

  vsd <- vst(dds, blind=FALSE)
  pca_data <- plotPCA(vsd, intgroup=c('condition'), returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  timepoints_used<-my_samp_dta$timepoint[order(match(my_samp_dta$sample,pca_data$name))]
  pca_data$timepoint<-factor(timepoints_used,levels=unique(timepoints_used))

  pca_data$group<-factor(pca_data$group,levels=c(group_order))


  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group, shape=timepoint,label = name)) +
    geom_label_repel(aes(PC1, PC2, label = name), fontface = 'bold',
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.5, "lines"),
                     segment.color = 'grey50') +
    scale_shape_manual(values=1:8) +
    geom_point(size=2, stroke = 2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))

  return(pca_plot)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' plotting function to plot the trajectories of genes in different groups
#'
#' @param ts_df dataframe containing the scaled value of individual genes for
#' every timepoint. dataframe is produced by\code{cluster_traj_merged_groups}
#' @param mean_df dataframe containing the mean of the group of genes for each
#' group being plotted. dataframe produced by \code{plot_merged_cluster_traj}
#' @param col_vect A vector containing colors whose names match the group names
#' in te ts_df and mean_df.
#'
#' @return The ggplot2 object for the plot
#'
#' @export
plot_merged_cluster_traj<-function(ts_df,mean_df,col_vect){
  plt <- ggplot(ts_df, aes(y = trans_mean , x = timepoint, color = group))+
    scale_color_manual(values=col_vect) +
    geom_line(aes(group = gene_id), alpha = 0.1) +
    geom_point() +
    geom_line(
      data = mean_df, lwd = 2, color = 'grey50',
      aes(group = group)
    ) +
    facet_wrap(~labels, scales = 'free_x', ncol = 1) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

  plt <- plt + theme(strip.text.x = element_text(size = 35))
  plt <- plt + theme(axis.text.x = element_text(size=25))
  plt <- plt + theme(axis.text.y = element_text(size=25))
  return(plt)
}
