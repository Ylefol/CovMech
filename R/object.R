# library(DESeq2)
# library(limma)
packages_for_loading<-c('DESeq2','limma')
suppressPackageStartupMessages(lapply(packages_for_loading, require, character.only = TRUE))


# Class creation  -------------------


#' The TimeSeries class
#'
#' The TimeSeries class is the main class for the timeseries pipeline
#' It is used to store the input data along with the final processed data and
#' any intermediate results that may be used for plot creation
#'
#' @slot count_matrix A list used to store the raw and normalized count matrices
#' @slot sample_data A dataframe created by the user indicating the grouping, timepoints, 
#' and replicates to which the samples belong 
#' @slot group_names A vector for the group names used, should be in order for 
#' which the differential expression comparisons will occur, where the first 
#' should be the experiment and second should be the control
#' @slot group_colors A named vector indicating which colors should be associated to which group
#' @slot DE_method Either limma or DESeq2 to indicate which differential expression method is used
#' @slot DE_p_filter Character of paj or pvalue to use in the establishment of 
#' significant differential expressed genes
#' @slot DE_p_thresh A numeric value used along with the pvalue filter to establish
#' significance. Genes will be considered significant if they are below this threshold
#' @slot DE_l2fc_thresh A numeric value for the log2FolChange threshold to establish 
#' significance where the absolute value of the genes must be greater than the threshold
#' @slot Gpro_org The organism used in gprofiler firendly format
#' @slot DESeq2_obj The normalized DESeq2 object
#' @slot limma_object The EList limma object
#' @slot DE_results A list of results for the different differential expression experiments
#' performed
#' @slot PART_L2FC_thresh A integer indicating the log(2)foldchange threshold for genes
#' to be PART clustered
#' @slot PART_results A list of results for the PART clustering analysis
#' @slot sem_sim_org A string indicating the annotation DBI organism to use
#' @slot Gprofiler_results A list of Gprofiler results for the PART clusters
#'
#' @name TimeSeries_Object-class
#' @rdname TimeSeries_Object-class
#' @concept objects
#' @exportClass TimeSeries_Object
#'
TimeSeries_Object<-setClass(
  Class='TimeSeries_Object',slots=list(
    count_matrix='list',
    sample_data='data.frame',
    group_names='vector',
    group_colors='vector',
    DE_method='character',
    DE_p_filter='character',
    DE_p_thresh='numeric',
    DE_l2fc_thresh='numeric',
    Gpro_org='character',
    DESeq2_obj='DESeqDataSet',
    limma_object='EList',
    DE_results='list',
    PART_l2fc_thresh='numeric',
    PART_results='list',
    sem_sim_org='character',
    Gprofiler_results='list'))


# Class related functions  -------------------
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A simple function to read and subset the sample_data created by the user
#' 
#'
#' @param path The path that will be read to obtain the sample data
#' @param group_names The group names selected by the user
#'
#' @return the subsetted sample file based on the selected groups
#'
#' @examples
#' 
#'
#' @export
prep_sample_data<-function(path, group_names){
  sample_file<-read.csv(path)
  sample_file<-sample_file[sample_file$group %in% group_names,]
  
  return(sample_file)
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which prepares and ads the raw count matrix to the timeseries object
#' 
#' The function also differs in the preparation method depending onthe type of
#' differential gene expression method selected (DESeq2 or limma) as this affects
#' the input data. Where DESeq2 (RNAseq data) is stored as individual count file while
#' for limma the function expects a csv file with all samples and genes.
#' 
#'
#' @param time_object A timeseries object
#' @param path_to_data The path towards the count files.
#' @param limma_from_raw Indicates if raw data is provided for microarrays or if an 
#' Elist was given
#' @param limma_id_replace The name of gene id's to use from the 'genes' dataframe
#' in the Elist
#'
#' @return The timeseries object with the raw count matrix added to it
#'
#' @examples
#' 
#'
#' @export
create_raw_count_matrix<-function(time_object,path_to_data=NULL,limma_id_replace='GeneName'){
  
  groups<-time_object@group_names
  #Ensures that the order will follow the grouping order
  selected_samples_1<-time_object@sample_data$sample[time_object@sample_data$group %in% groups[1]]
  selected_samples_2<-time_object@sample_data$sample[time_object@sample_data$group %in% groups[2]]
  selected_samples<-c(selected_samples_1,selected_samples_2)
  
  #Prepare the matrix according to the differential expression method (affects input)
  if (time_object@DE_method=='DESeq2'){
    final_counts<-prep_RNAseq_matrix(path_to_data,selected_samples)
    slot_storage<-'raw'
  }else if (time_object@DE_method=='limma'){
    if(endsWith(path_to_data,'.rds')){
      time_object@limma_object<-readRDS(my_path_data)
    }else{
      time_object@limma_object<-process_microarr_dta_limma(my_path_data)
    }
    
    final_counts<-prep_limma_matrix(Elist_obj=time_object@limma_object,replace_rows_with = limma_id_replace)
    slot_storage<-'norm' #Stored as norm as the data is expected to be normalized already
  }
  #Re-organize samples to be in order of groups
  final_counts<-final_counts[,selected_samples]
  
  
  time_object@count_matrix[[slot_storage]]<-final_counts
  return(time_object)
}

#' Function to read and process microarray data using the limma package
#' 
#' The function reads in the microarray data, removes the control probes, performs
#' background correction, normalization, and averages over irregular replicate probes (avereps)
#' 
#' @param micro_arr_path String indicating the path to a folder containing all and only 
#' microrna text
#' @param micro_arr_source The source of the microarray data. Choices can be seen using
#' limma documentation for the function \code{read.maimages}
#' @param green.only Boolean for the green only parameter of \code{read.maimages}
#' @param backg_corr_meth The background correction method with \code{backgroundCorrect}
#' @param back_corr_offset The offset integer to use with \code{backgroundCorrect}
#' @param norm_arrays_meth The normalization method to use with \code{normalizeBetweenArrays}
#' @param ID_used The gene IDs used for the Elist. The name given must match a column in the 
#' 'genes' dataframe of the data. If set to NULL, the first column is taken
#' 
#' @return The processed and normalized Elist object
#' 
process_microarr_dta_limma<-function(micro_arr_path,micro_arr_source='agilent',green.only=TRUE,
                                     backg_corr_meth='normexp',back_corr_offset=16,
                                     norm_arrays_meth='quantile',ID_used='GeneName'){
  #Put file names in vector
  files<-list.files(micro_arr_path)
  raw_files_path <- paste(micro_arr_path, files, sep='/') # add path to filenames
  
  # load data (see section 4.5 in the limma user guide).
  rawdata <- read.maimages(file=raw_files_path, source=micro_arr_source, green.only=green.only)
  
  #Rename columns to not include the path to the files
  new_colnames<-gsub(x = colnames(rawdata$E),pattern = paste0(micro_arr_path,'/'),replacement = '')
  colnames(rawdata$E)=new_colnames
  
  #Remove control
  ctrl_idx <- rawdata@.Data[[4]]$ControlType != 0
  rawdata_noctrl <- rawdata[!ctrl_idx,]
  
  #Background correction and normalization
  data_bckcor <- backgroundCorrect(rawdata_noctrl, method=backg_corr_meth, offset=back_corr_offset)
  norm_data <- normalizeBetweenArrays(data_bckcor, method=norm_arrays_meth) 
  
  #Create gene ID vector and average the reps
  if(is.null(ID_used)==T){
    ID_vect<-norm_data$genes[,1]
  }else{
    ID_vect<-norm_data$genes[[ID_used]]
  }
  final_dta <- avereps(norm_data, ID=ID_vect)
  
  return(final_dta)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which extracts and modifies (if needed) a data matrix from a Elist object
#' 
#' The function can replace ronames of the matrix by one contained in the 'gene'
#' container of the Elist
#' 
#'
#' @param Elist_obj A Elist object
#' @param replace_rows_with Character indicating with column (if any) will be used
#' to replace the rownames of the matrix in the Elist
#'
#' @return count_matrix 
#'
#' @examples
#' 
#'
#' @export
prep_limma_matrix<-function(Elist_obj,replace_rows_with=NULL){

  count_matrix=Elist_obj$E
  #Replaces rownames with a selection from the genes dataframe of the E list
  if (is.null(replace_rows_with)==F){
    row.names(count_matrix)=Elist_obj$genes[[replace_rows_with]]
  }  
  
  return(count_matrix)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which prepares a count matrix form RNAseq count files
#' 
#' Function requires the path to the individual count files and the samples to 
#' be included in the matrix. Files will be read based on the sample they represent
#' the values from the different files are merged into a matrix and returned
#' 
#' 
#'
#' @param path_to_counts Path to the individual RNAseq count files
#' @param selected_samples The samples to be read
#'
#' @return final_counts The formatted count matrix
#'
#' @examples
#' 
#'
#' @export
prep_RNAseq_matrix<-function(path_to_counts,selected_samples){
  final_counts<-data.frame(NULL)
  for(file in list.files(path_to_counts)){
    sample_name<-strsplit(file,'\\.')[[1]][1]
    if (sample_name %in% selected_samples){
      if(endsWith(file,'.csv')==T){
        temp_df<-read.csv(paste0(path_to_counts,'/',file),header=T)
      }else{
        temp_df<-read.delim(paste0(path_to_counts,'/',file),header=F)
      }
      colnames(temp_df)=c('gene_id',sample_name)
      if (nrow(final_counts)==0){
        final_counts<-temp_df
      }else{
        final_counts<-merge(final_counts,temp_df,by='gene_id')
      }
    }
  }
  #Remove any non-genes (starts with underscore)
  final_counts<-final_counts[!startsWith(final_counts$gene_id,'_'),]
  row.names(final_counts)=final_counts$gene_id
  final_counts<-final_counts[,colnames(final_counts)!='gene_id']
  
  return(final_counts)
}




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which loads example data from a alternate time series pipeline in order
#' to compare the results. The function loads and formats this data to the specifications
#' of this pipeline.
#' 
#' @param p_val_filter_type Character of paj or pvalue to use in the establishment of 
#' significant differential expressed genes
#' @param p_thresh A numeric value used along with the pvalue filter to establish
#' significance. Genes will be considered significant if they are below this threshold
#' @param l2fc_thresh A numeric value for the log2FolChange threshold to establish 
#' significance where the absolute value of the genes must be greater than the threshold
load_cop1_example_data<-function(p_val_filter_type='padj',p_thresh=0.05,l2fc_thresh=1){
  urls <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/",
    "?acc=GSE114762&format=file&file=",
    c(
      "GSE114762_raw_counts.csv.gz",
      "GSE114762_gene_data.csv.gz",
      "GSE114762_sample_info.csv.gz"
    )
  )
  
  library(BiocFileCache)
  library(tibble)
  library(dplyr)
  library(tidyverse)
  bfc <- BiocFileCache(ask = FALSE)
  
  cnts <- read_csv(bfcrpath(rnames = urls[1])) %>%
    remove_rownames() %>%
    as.matrix()
  rownames(cnts) <- cnts[,1]
  cnts <- cnts[,-1]
  
  gene.data <- read_csv(bfcrpath(rnames = urls[2])) %>%
    as.data.frame() %>%
    remove_rownames() 
  rownames(gene.data) <- gene.data[,1]
  gene.data <- gene.data[,-1]
  
  # row.names(cnts)=gene.data$symbol
  # row.names(cnts)=make.unique(gene.data$symbol)
  library("org.Mm.eg.db")
  my_org_sem_sim='org.Mm.eg.db'
  gene_id<-mapIds(org.Mm.eg.db,keys=row.names(cnts), column="SYMBOL",keytype="ENTREZID")
  annot_df<-data.frame(gene_id)
  annot_df$gene_id <- ifelse(is.na(annot_df$gene_id), row.names(annot_df), annot_df$gene_id)
  row.names(cnts)=annot_df$gene_id
  
  pheno.data <- read_csv(bfcrpath(rnames = urls[3])) %>%
    as.data.frame() %>%
    remove_rownames() 
  rownames(pheno.data) <- pheno.data[,1]
  pheno.data <- pheno.data[,-1]
  pheno.data<-pheno.data[,c('sample','group','replicate','time')]
  colnames(pheno.data)=c('sample','group','replicate','timepoint')
  
  graph_vect<-c("#e31a1c","#1f78b4")
  names(graph_vect)=c('Loxp','WT')
  time_object <- new('TimeSeries_Object',sample_data=pheno.data,
                   group_names=c('Loxp','WT'),group_colors=graph_vect,DE_method='DESeq2',
                   DE_p_filter=p_val_filter_type,DE_p_thresh=p_thresh,DE_l2fc_thresh=l2fc_thresh,
                   Gpro_org='mmusculus')
  
  time_object@count_matrix[['raw']]<-cnts
  
  return(time_object)
}


#' #' The custom_DE_results class
#' #'
#' #' Class used to store differential expression results, split into conditional
#' #' and temporal. Where each can store any amount of differential expression expriments
#' #' performed. Each experiment will contain the raw DE object (a deseq2 object for 
#' #' RNAseq and a MArrayLM object containing the eBayes results). They will also contain
#' #' the raw and significant differential expression data. The significance parameters
#' #' are defined by the parameters in the time series object
#' #'
#' #' @slot conditional results for conditional differential gene expression
#' #' @slot temporal results for temporal differential gene expression
#' #'
#' #' @name custom_DE_results-class
#' #' @rdname custom_DE_results-class
#' #' @concept objects
#' #' @exportClass custom_DE_results
#' #'
#' custom_DE_results<-setClass(
#'   Class='custom_DE_results',slots=list(
#'     conditional='list',
#'     temporal='list'))
#' 
#' 
#' 
#' #' The PART_results class
#' #' 
#' #' Class which contains the various elements of the PART clustering analysis
#' #'
#' #' @slot part_matrix The matrix subset used in the PART analysis, values are transformed
#' #' @slot cluster_info The hclust and calculated clusters obtained by PART clustering
#' #' @slot part_data A data frame containing the genes from the PART matrix
#' #' @slot cluster_map A dataframe containing the list of genes, associated cluster
#' #' and associated hex code for color representation.
#' #'
#' #' @name PART_results-class
#' #' @rdname PART_results-class
#' #' @concept objects
#' #' @exportClass PART_results
#' #'
#' PART_results<-setClass(
#'   Class='PART_results',slots=list(
#'     part_matrix='matrix',
#'     cluster_info='list',
#'     part_data='data.frame',
#'     cluster_map='data.frame'))


