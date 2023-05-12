library(clusterGenomics)
library(DESeq2)
library(viridis)
library(Biobase)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(TimeSeriesExperiment)
library(UpSetR)
library(limma)
library(ComplexHeatmap)


# data retrieval functions  -------------------

#' Function which loads the dataset from the TimeSeriesExperiment tutorial
#'
#' The dataset represents cells from six mice. Cells are divided in two equal groups, 
#' in one group the Cop1 gene is KO'd using tamoxifen. All subjects were then subject 
#' to LPS treatment to induce an inflammatory response. 6 time points were bulk 
#' RNAseq'd, time 0 which is before LPS, then the time points of 2.5, 4, 6, 9, 
#' and 13 hours after addition of LPS
#' 
#' @return TS_object A time series object containing the example dataset
#'
#' @examples
#' TS_dat <- load_example_data()
#'
#' @export
load_example_data <-function(){
  #Load the urls for the datasets
  urls <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/",
    "?acc=GSE114762&format=file&file=",
    c(
      "GSE114762_raw_counts.csv.gz",
      "GSE114762_gene_data.csv.gz",
      "GSE114762_sample_info.csv.gz"
    )
  )
  
  
  #Load library and download dataset
  library(BiocFileCache)
  bfc <- BiocFileCache(ask = FALSE)
  
  cnts <- read_csv(bfcrpath(rnames = urls[1])) %>%
    remove_rownames() %>%
    column_to_rownames("X1") %>%
    as.matrix()
  
  gene.data <- read_csv(bfcrpath(rnames = urls[2])) %>%
    as.data.frame() %>%
    remove_rownames() %>%
    column_to_rownames("X1")
  
  gene.data <- gene.data[1:2]
  print(head(gene.data))
  pheno.data <- read_csv(bfcrpath(rnames = urls[3])) %>%
    as.data.frame() %>%
    remove_rownames() %>%
    column_to_rownames("X1")
  

  TS_object <- TimeSeriesExperiment(
    assays = list(cnts), # seems to be count data
    rowData = gene.data, # gtf-like information (feature, symbol, size, type, description)
    colData = pheno.data, # phenotypic data (sample, group, individual,replicate,time,treatment,label)
    timepoint = "time",
    replicate = "replicate",
    group = "group"
  )
  
  rowData(TS_object)$feature=rowData(TS_object)$symbol
  return (TS_object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which loads the dataset from the LEMONAID project
#'
#' Retrieves data from a pre-determined location based on the inputted groups
#' 
#' @return TS_object A time series object containing the example dataset
#'
#' @examples
#'
#' @export
TS_obj_LEMONAID_data <- function(group_select){
  setwd("~/A_Projects/EpiGen/R_Work_Folder/rnaseq_pipeline_yl")
  
  
  count_list <- list.files('side_work/AID_TS_folder/raw_counts_TS/')
  
  for (file in count_list){
    my_file <- read.table(paste0('side_work/AID_TS_folder/raw_counts_TS/',file),sep='\t',header = F)
    file_name<-strsplit(file,'.counts')[[1]]
    file_name<-gsub(x = file_name,pattern = '-',replacement = '')
    colnames(my_file)=c('gene_id',file_name)
    if (exists('final_counts')==F){
      final_counts <- my_file
    }else{
      final_counts <- merge(final_counts,my_file,by='gene_id')
    }
  }
  
  rownames(final_counts)=final_counts$gene_id
  final_counts <-final_counts[,-1]
  
  
  
  #Remove first five rows, it relates to non aligned counts
  final_counts <- tail(final_counts, -5)

  count_names <- colnames(final_counts)
  my_counts <- final_counts
  my_pheno <- read.table('side_work/AID_TS_folder/id_explanation_for_R_no_27_70.csv',sep=',',header=T,row.names = 1)
  group_select<-c('Control',group_select)
  my_pheno<-my_pheno[my_pheno$group%in%group_select,]
  
  my_counts<-my_counts[,my_pheno$sample]
  symbol <- as.vector(rownames(final_counts))
  
  library(AnnotationDbi)
  library("org.Hs.eg.db")
  selected_org_obj=org.Hs.eg.db
  feature <- as.vector(mapIds(selected_org_obj, keys = symbol, keytype = "SYMBOL", column="ENTREZID"))
  
  my_gene_info <- data.frame(symbol, feature,row.names = symbol)
  #re-order
  col_order <- as.vector(rownames(my_pheno))
  my_counts <- my_counts[, col_order]

  TS_object <- TimeSeriesExperiment(
    assays = list(my_counts), # seems to be count data
    rowData = my_gene_info, # gtf-like information (feature, symbol, size, type, description)
    colData = my_pheno, # phenotypic data (sample, group, individual,replicate,time,treatment,label)
    timepoint = "time",
    replicate = "replicate",
    group = "group"
  )
  TS_object@elementMetadata@listData[["feature"]] <- my_gene_info$feature
  return(TS_object)
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which loads the Robin's AID murine dataset
#'
#' Retrieves data from a pre-determined location based on the inputted groups
#' 
#' @return TS_object A time series object containing the example dataset
#'
#' @examples
#'
#' @export

TS_obj_Robins_AID_data <- function(){
  count_list <- list.files('side_work/Robin_AID/merged_counts/')
  
  for (file in count_list){
    my_file <- read.table(paste0('side_work/Robin_AID/merged_counts/',file),sep=',',header = T, row.names = 1)
    if (exists('final_counts')==F){
      final_counts <- my_file
    }else{
      final_counts <- merge(final_counts,my_file,by='gene_id')
    }
  }
  rownames(final_counts)=final_counts$gene_id
  final_counts <-final_counts[,-1]
  
  count_names <- colnames(final_counts)
  new_names <-c()
  for (name in count_names){
    reduced_name <- strsplit(name,'_')
    new_names <- c(new_names,reduced_name[[1]][1])
  }
  
  colnames(final_counts)=new_names
  
  #Duplicate each experiment to create 'replicates'
  for (exp in colnames(final_counts)){
      final_counts[paste0(exp,'_2')] <- final_counts[exp]
    
    # if (endsWith(x = exp,suffix = 'C')){
    #   # print(exp)
    #   final_counts[paste0(exp,'_2')] <- final_counts[exp]
    # }
  }
  
  
  my_counts <- final_counts
  my_pheno <- read.table('side_work/Robin_AID/sample_information_revamp_2.csv',sep=',',header=T,row.names = 1)
  
  symbol <- as.vector(rownames(final_counts))
  View(my_counts)
  library(AnnotationDbi)
  library("org.Mm.eg.db")
  selected_org_obj=org.Mm.eg.db
  feature <- as.vector(mapIds(selected_org_obj, keys = symbol, keytype = "SYMBOL", column="ENTREZID"))
  
  my_gene_info <- data.frame(symbol, feature,row.names = symbol)
  #re-order
  col_order <- as.vector(rownames(my_pheno))
  my_counts <- my_counts[, col_order]
  
  TS_object <- TimeSeriesExperiment(
    assays = list(my_counts), # seems to be count data
    rowData = my_gene_info, # gtf-like information (feature, symbol, size, type, description)
    colData = my_pheno, # phenotypic data (sample, group, individual,replicate,time,treatment,label)
    timepoint = "time",
    replicate = "replicate",
    group = "group"
  )
  # TS_object@elementMetadata@listData[["feature"]] <- my_gene_info$feature
  return(TS_object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which loads the COVID data
#'
#' Retrieves data from a pre-determined location based on the inputted groups
#' @param use_L2FC Wether or not to use the raw micro-array data or formatted data 
#' @return TS_object A time series object containing the example dataset
#'
#' @examples
#'
#' @export
load_covid_data <- function(use_L2FC=F){
  if (use_L2FC==T){
    final_counts <- read.csv('side_work/covid/formatted_data/formatted_TS_with_healthy_covid_data.csv',check.names = FALSE)
  }else{
    final_counts <- read.csv('side_work/covid/formatted_data/formatted_TS_raw.csv',check.names = FALSE)
  }
  
  rownames(final_counts)=final_counts$GeneName
  gene_symbols <- as.vector(final_counts$GeneName)
  feature_vect <- as.vector(final_counts$EntrezID)
  
  my_gene_info <- data.frame(gene_symbols, feature_vect,row.names = gene_symbols)
  
  gene_info_conv <- my_gene_info[is.na(my_gene_info$feature_vect)==T,]
  gene_info_conv <- gene_info_conv[startsWith(gene_info_conv$gene_symbols,'ENST'),]
  gene_info_conv<- gene_info_conv %>% 
    rename(
      transcript_id = gene_symbols
    )
  #Load gtf
  gtf <- rtracklayer::import('data_files/gencode.v35.annotation.gtf')
  gtf=as.data.frame(gtf)
  #Find ENSTs while removing decimals
  sub_gtf <- gtf[gsub('\\..','',gtf$transcript_id) %in% gene_info_conv$transcript_id,]
  sub_gtf <- sub_gtf[!duplicated(sub_gtf$gene_name),]
  sub_gtf <- sub_gtf[,c('gene_name','transcript_id')]
  #Rename to be able to merge and identify new symbols
  sub_gtf<- sub_gtf %>% 
    rename(
      new_gene_symbols = gene_name
    )
  #Removal of decimals
  sub_gtf$transcript_id <- gsub('\\..','',sub_gtf$transcript_id)
  sub_gtf$new_gene_symbols <- gsub('\\..','',sub_gtf$new_gene_symbols)
  merged_gtf_conv <- merge(gene_info_conv,sub_gtf,by='transcript_id')

  #Find the features using the new gene symbols
  library(AnnotationDbi)
  library("org.Hs.eg.db")
  selected_org_obj=org.Hs.eg.db
  merged_gtf_conv$new_features <- as.vector(mapIds(selected_org_obj, keys = merged_gtf_conv$new_gene_symbols, keytype = "SYMBOL", column="ENTREZID"))

  #Overwrite the necessary columns of the count matrix, recreate my_gene_info with
  for (i in 1:length(gene_symbols)){
    if (gene_symbols[i] %in% merged_gtf_conv$transcript_id){
      feature_vect[i] <- merged_gtf_conv$new_features[merged_gtf_conv$transcript_id==gene_symbols[i]]
      gene_symbols[i] <- merged_gtf_conv$new_gene_symbol[merged_gtf_conv$transcript_id==gene_symbols[i]]
    }
  }
  gene_symbols <- make.unique(gene_symbols,sep='.')
  my_gene_info <- data.frame(gene_symbols, feature_vect,row.names = gene_symbols)
  
  #Replace relevant columns in count matrix
  final_counts$GeneName <- gene_symbols
  rownames(final_counts)=final_counts$GeneName

  final_counts <- final_counts[,c(-1,-2)]
  
  my_counts <- final_counts
  # my_pheno <- read.csv('side_work/covid/formatted_data/formatted_covid_clin_data_trips.csv',check.names=FALSE,row.names = 1)
  my_pheno <- read.csv('side_work/covid/formatted_data/formatted_covid_clin_data_trips_M_60_below.csv',check.names=FALSE,row.names = 1)
  # my_pheno <- read.csv('side_work/covid/formatted_data/clin_trips_with_cluster_grouping.csv',check.names=FALSE,row.names = 1)
  # my_pheno <- read.csv('side_work/covid/formatted_data/clin_trips_with_cluster_grouping_pure.csv',check.names=FALSE,row.names = 1)

  my_counts <- my_counts[colnames(my_counts) %in% rownames(my_pheno)]
  
  #re-order
  col_order <- as.vector(rownames(my_pheno))

  my_counts <- my_counts[, col_order]
  
  my_counts <- round(my_counts)
  
  TS_object <- TimeSeriesExperiment(
    assays = list(my_counts), # seems to be count data
    rowData = my_gene_info, # gtf-like information (feature, symbol, size, type, description)
    colData = my_pheno, # phenotypic data (sample, group, individual,replicate,time,treatment,label)
    timepoint = "time",
    replicate = "replicate",
    group = "group"
  )
  TS_object@elementMetadata@listData[["symbol"]] <- my_gene_info$gene_symbols
  TS_object@elementMetadata@listData[["feature"]] <- my_gene_info$feature_vect
  return(TS_object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which uses limma data to create a time series object
#'
#' Loads the limma data from a given path. The objects in the path are expected
#' to be called 'rna_biop_dataa' and be an E-list object.
#' sample data is also provided via a path and is expected to contain a column called
#' 'sample' with values/names matching the sample names within the E-list
#' 
#' This serves as a means to select which samples from the Elist are to be used in
#' the time series object
#' 
#' @param path_name_Elist path to the Elist object/Rdata
#' @param sample_path_name path to the csv file containing a column of samples
#' 
#' @return TS_object A time series object containing the example dataset
#'
#' @examples
#'
#' @export
limma_to_TS <-function(path_name_Elist,sample_path_name){
  load(path_name_Elist)
  count_matrix<-rna_biop_dataa$E
  TS_sample_data<-read.csv(sample_path_name,row.names = 1)
  count_matrix<-count_matrix[,TS_sample_data$sample]
  row.names(count_matrix)=rna_biop_dataa$genes$GeneName
  
  gene_info<-rna_biop_dataa$genes
  gene_info<-gene_info[,c('GeneName','EntrezID')]
  colnames(gene_info)=c('symbol','feature')
  gene_info$symbol<-make.unique(gene_info$symbol)
  row.names(gene_info)=gene_info$symbol
  
  TS_object <- TimeSeriesExperiment(
    assays = list(count_matrix),
    rowData = gene_info,
    colData = TS_sample_data,
    timepoint = "time",
    replicate = "replicate",
    group = "group"
  )
  TS_object@elementMetadata@listData[["feature"]] <- gene_info$feature
  return(TS_object)
}


# pipeline functions  -------------------


#' A wrapper function which formats the normalizes the count data from a time series 
#' experiment using DESeq2's normalization
#' 
#' The function takes the timeseries object along with the groups it contains,
#' it then splits the two groups into two separate files which are inputted in the
#' DESeq2 normalization function, another function is called to store the normalized
#' counts into the time object itself, the time object is then returned
#'
#' 
#' @param time_object A time series object
#' @param groups The name of the groups in the TS object (must be exact copies)
#' 
#' @param time_object The time object with the normalized count matrix included
#'
#' @examples
#' TS_dat <- load_example_data()
#' TS_object <- normalize_timeSeries_with_deseq2(time_object=TS_object,groups=c('WT','Loxp'))
#'
#' @export
normalize_timeSeries_with_deseq2 <- function(time_object, groups=c('control','experiment')){
  for (item in groups){
    file_sample <- colData(time_object)[colData(time_object)$group==item,]
    file_sample <- file_sample$sample
    file <- assays(time_object)$raw[,colnames(assays(time_object)$raw) %in% file_sample]
    file <-as.data.frame(file)
    file<-add_column(file, Genes = rownames(assays(time_object)$raw), .before = 1)
    
    if (exists('file_1')==T){
      file_2 <-file
    }else{
      file_1<-file
    }
  }
  norm_counts <- deseq2_normalization(file_1,file_2,groups[1],groups[2])
  #Re-orders norm as to have the same order as the raw assay, prevents re-ordering issues
  #Downstream in the pipeline
  norm_counts<-norm_counts[,colnames(assays(time_object)$raw)]
  time_object <- normalizeData_mod(time_object,norm_counts)
  return(time_object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A function which normalizes count files as per DESeq2's normalization
#' 
#' The function takes in two files and two categories/conditions
#' The files are filtered and merged to then be transformed into a deseq2 object
#' The deseq2 object is normalized and returned in a count matrix format.
#'
#' 
#' @param File_1 The dataframe file for category/condition 1, must contain a 'Genes' column
#' @param File_2 The dataframe file for category/condition 2, must contain a 'Genes' column
#' @param File_1_category The name of the first condition
#' @param File_2_category The name of the second condition
#' 
#' @return my_counts The normalized count matrix
#'
#' @examples
#' norm_counts <- deseq2_normalization(file_1,file_2,groups[1],groups[2])
#' 
#' @export
deseq2_normalization <- function(File_1, File_2, File_1_category, File_2_category){
  
  # preprocess both files, removal of NAs and counts with values of 0.0
  File_1<-File_1[complete.cases(File_1),]
  # File_1<-File_1[!(File_1$Genes=="0.0"),]
  
  File_2<-File_2[complete.cases(File_2),]
  # File_2<-File_2[!(File_2$Genes=="0.0"),]
  
  #Counts the number of columns/patients in each file
  num_1 <- ncol(File_1)-1
  num_2 <- ncol(File_2)-1
  
  #Remove first column of file that is being added
  #These are the counts that are will be merged to another
  new_file_2<- File_2[,-1,drop=FALSE]
  nams=File_1[['Genes']]
  
  #Create 'merged' files
  merged_files <- cbind(File_1,new_file_2)
  merged_files <- merged_files[,-1]
  
  #converts it to numeric
  merged_files <- data.frame(sapply(merged_files, as.numeric),check.names=F, row.names = rownames(merged_files))
  
  #Creates the matrices
  merged_files <- as.matrix(merged_files)
  
  #Assign Condition to all of the columns/samples/patients
  condition <- factor(c(rep(File_1_category, num_1),rep(File_2_category,num_2)))
  
  #Create a coldata frame and instantiate the DESeqDataSet
  col_data <- data.frame(row.names=colnames(merged_files),condition)
  dds <- DESeqDataSetFromMatrix(countData=merged_files, colData=col_data, design=~condition)
  dds = estimateSizeFactors(dds)
  
  my_counts <- as.data.frame(counts(dds,normalized=TRUE))
  my_counts<- na.omit(my_counts)
  
  return(my_counts)
  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A function which adds the normalized count data to the TS_object
#' 
#' @param object The time series object
#' @param the_counts The The normalized count matrix
#' 
#' @return object The time series object with the added normalized count matrix
#'
#' @examples
#'norm_counts <- deseq2_normalization(file_1,file_2,groups[1],groups[2])
#'time_object <- normalizeData_mod(time_object,norm_counts)
#' 
#' @export
normalizeData_mod <- function(object, the_counts) 
{
  curr.assays <- assays(object)
  raw.data <- normalized.data <- curr.assays$raw
  normalized.data <- as.matrix(the_counts)
  if (is(object, "TimeSeriesExperiment") ) {
    curr.assays$norm <- normalized.data
    slot(object, name = "assays", check = TRUE) <- Assays(curr.assays)
  }
  return(object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' A wrapper function which filters a time series data based on minimum mean gene expression
#' 
#' The function takes in the time series object and an integer value.
#' The integer is used to set the minimum mean counts for a gene in any of the groups
#' The filter is done on the normalized data, the mean is in relation to the means
#' of all the replicates put together for a single gene at a single timepoint.
#' 
#' @param object A time series object containing normalized count data
#' @param filter_value The integer value for the minimum mean filter
#'
#' 
#' @return object A minimum mean fitlered time series object
#'
#' @examples
#' TS_object <-minimum_mean_TS_filter(TS_object,5)
#'
#' @export
minimum_mean_TS_filter <- function (object, filter_value=5){
  min_mean_cpm <- filter_value
  group_cpm_means <- data.frame(row.names = rownames(object))
  norm.cnts <- assays(object)$norm
  for(g in unique(TimeSeriesExperiment::groups(object))) {
    g_cnts <- norm.cnts[ , which(TimeSeriesExperiment::groups(object) == g)]
    group_cpm_means[, g] <- apply(g_cnts, 1, mean)
  }
  group_cpm_max <- apply(as.vector(group_cpm_means), 1, max)
  genes_expressed <- rownames(object)[group_cpm_max > min_mean_cpm]
  object <- filterFeatures(object, genes_expressed)
  return(object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Creates a cluster map based on csv data from a heatmap
#' 
#' Function reads a csv from a heatmap plot, it then formats the csv information
#' into a cluster map which can be used to plot specific cluster time series plots
#' for individual genes.
#' 
#' The function assumes that genes are shown as gene symbol in the csv and converts
#' them to feature as it is necessary for the clsuter time series plots.
#' 
#' @param csv_path The path to the csv that will be used to create the cluster map
#'
#' 
#' @return cluster_map A cluster map
#'
#' @examples
#' plotHeatmap_part_mod(TS_object, num.feat = 300, save_name = 'TS_results/top_300_variable_genes.csv')
#' cluster_map <- create_cluster_map_from_csv(TS_object,csv_path='TS_results/top_300_variable_genes.csv')
#'
#' @export
create_cluster_map_from_csv <- function(csv_path){
  top_300<-read.csv(csv_path)
  cluster_map <- top_300[,c(1,2,3)]
  cluster_map$gene_cluster <- paste("C", cluster_map$gene_cluster,sep='')
  colnames(cluster_map)=c('feature','cluster','cluster_col')
  cluster_map$feature <- as.character(cluster_map$feature)
  return(cluster_map)
}


# Differential Expression functions -------------------


#' Function which prepares a timepoint model matrix for limma datasets
#' 
#' The function creates either a conditional matrix (do_temporal==FALSE) where
#' the model matrix will compare both conditions ex: WT vs Knockout
#' Or it will create a model for a temporal comparison, where it will compare the
#' two timpoints contained in 'groups_in_ts'
#' 
#' The function is meant to be used within a wrapper function
#' 
#' @param time_object A time series object
#' @param target_tp Vector containing the two timepoints to compare
#' @param groups_in_ts Vector of the two groups/conditionals
#' @param microarr_dta An E list object
#' @param do_temporal boolean indication if the model matrix should be temporal or conditional
#'
#' returns list(mm,subset_samples) A list contianing the model matrix and the samples
#' used to create the model matrix
#'
#' @examples
#' 
#'
#' @export
prep_tp_matrix <- function(time_object,target_tp,groups_in_ts,microarr_dta,do_temporal=FALSE){
  if(do_temporal==FALSE){
    group_ctrl<-colData(time_object)[colData(time_object)$group==groups_in_ts[1],]$sample
    group_exp<-colData(time_object)[colData(time_object)$group==groups_in_ts[2],]$sample
    tp_samples<-names(which(timepoints(time_object)==target_tp))
    group_1_tp<-group_exp[group_exp %in% tp_samples]
    group_2_tp<-group_ctrl[group_ctrl %in% tp_samples]
  }else{
    #Account for order reversal (use higher timepoint as experiment)
    group_1_tp<-names(which(timepoints(time_object)==groups_in_ts[2]))
    group_2_tp<-names(which(timepoints(time_object)==groups_in_ts[1]))
  }


  subset_samples<-c(group_1_tp,group_2_tp)
  
  # tp_limma_subset <- rna_biop_dataa[,colnames(rna_biop_dataa$E) %in% subset_samples]
  g1_vect<-gsub(pattern = F,replacement = 0,x = colnames(microarr_dta$E) %in% group_1_tp)
  g1_vect<-gsub(pattern = T,replacement = 1,x = g1_vect)
  
  
  g2_vect<-gsub(pattern = F,replacement = 0,x = colnames(microarr_dta$E) %in% group_2_tp)
  g2_vect<-gsub(pattern = T,replacement = 1,x = g2_vect)
  matrix_df<-data.frame(Group1=g1_vect, Group2=g2_vect)
  mm <- model.matrix(~Group1+Group2,matrix_df)
  colnames(mm)=c('X.Intercept','G1','G2')
  
  return(list(mm,subset_samples))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Wrapper function to perform a DESeq2 differential expression analysis on a time
#' series object
#' 
#' The wrapper function can perform a conditional or temporal differential expression
#' analysis. The function uses DESeq2 and should therefore be used with RNAseq data
#' 
#' The function will create a folder within the TS_results folder and store the 
#' analysis results in a relevantly named folder
#' 
#' @param time_object A time series object
#' @param groups a vector containing the groups/conditionals
#' @param genes_of_interest A list with genes of interest, if any
#' @param species Code for the species being analyzed
#' @param do_temporal Boolean indicating if the analysis is temporal or conditional
#'
#'
#' @examples
#' 
#'
#' @export
wrapper_TS_DESeq2_DE <- function(time_object, groups=c('control','experiment'),genes_of_interest=c(),species='hg38',do_temporal=FALSE){

  
  if(do_temporal==FALSE){
    tp_iter<-length(unique(timepoints(time_object)))
    save_path='TS_results/DE_results_conditional'
  }else{
    tp_iter<-length(unique(timepoints(time_object)))-1
    save_path='TS_results/DE_results_temporal'
  }
  gene_num_list <- list()
  for (tp_idx in 1:tp_iter){
    
    if (do_temporal==FALSE){
      tp<-unique(timepoints(time_object))[tp_idx]
      cond_1<-groups[1]
      cond_2<-groups[2]
      exp_name<-paste0(cond_2,'_vs_',cond_1,'_TP_',tp)
      
      possible_samples_tp<-names(which(timepoints(time_object)==tp))
      samples_G1<-possible_samples_tp[possible_samples_tp %in% names(which(TimeSeriesExperiment::groups(time_object)==groups[1]))]
      samples_G2<-possible_samples_tp[possible_samples_tp %in% names(which(TimeSeriesExperiment::groups(time_object)==groups[2]))]
      
    }else{
      tp<-NULL
      my_timepoints<-unique(timepoints(time_object))
      cond_1=paste0('TP_',my_timepoints[tp_idx])
      cond_2=paste0('TP_',my_timepoints[tp_idx+1])
      exp_name<-paste0(cond_2,'_vs_',cond_1)
      
      samples_G1<-names(which(timepoints(time_object)==my_timepoints[tp_idx]))
      samples_G2<-names(which(timepoints(time_object)==my_timepoints[tp_idx+1]))

    }
    message("running: ", exp_name)
    tp_file_1<-as.data.frame(assays(time_object)$raw[,samples_G1])
    tp_file_1<-add_column(tp_file_1, Genes = rownames(assays(time_object)$raw), .before = 1)
    
    tp_file_2<-as.data.frame(assays(time_object)$raw[,samples_G2])
    tp_file_2<-add_column(tp_file_2, Genes = rownames(assays(time_object)$raw), .before = 1)


    gtf_p <- 'data_files/gencode.v35.annotation.gtf'
    #Groups are 'reversed' since control is group 1, and we want the results to be in regards to the experiment
    DE_info<-DESeq2_gtf_wrapper(tp_file_2,tp_file_1,cond_2,cond_1,gtf_path=gtf_p,
                                genes_of_interest=genes_of_interest,
                                result_location = save_path,time_point = tp,species=species)
    gene_num_list[[exp_name]] <- DE_info
  }
  df <- data.frame(matrix(unlist(gene_num_list), nrow=length(gene_num_list), byrow=T))
  rownames(df)=names(gene_num_list)
  colnames(df)=c('Number of Genes', "Number of Significant Genes")
  write.csv(df,paste0(save_path,"/DE_summary.csv"), row.names = TRUE)
  create_tables_genes_interest(genes_of_interest,names(genes_of_interest),pval_or_padj=p_val_filter_type,
                               path_to_DE = paste0(save_path,'/'),results_path = paste0(save_path,'/tables_genes_interest/'))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Wrapper function to perform a limma differential expression analysis on a time
#' series object
#' 
#' The wrapper function can perform a conditional or temporal differential expression
#' analysis. The function uses limma and should therefore be used with microarray data
#' 
#' The function will create a folder within the TS_results folder and store the 
#' analysis results in a relevantly named folder
#' 
#' @param time_object A time series object
#' @param Elist_path The path to the limma E list Rdata object
#' @param group_names a vector containing the groups/conditionals
#' @param genes_of_interest A list with genes of interest, if any
#' @param gtf_path The path to the gtf file if one is to be used
#' @param filter_choice The filter to use (padj or pvalue)
#' @param species Code for the species being analyzed
#' @param do_temporal Boolean indicating if the analysis is temporal or conditional
#'
#'
#' @examples
#' 
#'
#' @export
wrapper_limma_conditional_DE <-function(time_object,Elist_path,group_names,genes_of_interest=c(),
                                      gtf_path=NULL,filter_choice='padj',species='hg38',do_temporal=FALSE){
  load(Elist_path)
  if (do_temporal==FALSE){
    my_save_path='TS_results/DE_results_conditional/'
  }else{
    my_save_path='TS_results/DE_results_temporal/'
  }
  
  dir.create(my_save_path)
  dir.create(paste0(my_save_path,'Robjects_checkpoints'))
  
  
  #always assume 1 is the control (therefore will be put second in limma)
  if(do_temporal==FALSE){
    tp_iter<-length(unique(timepoints(time_object)))
  }else{
    tp_iter<-length(unique(timepoints(time_object)))-1
  }
  
  
  gene_num_list<-list()
  eb_res_list<-list()
  my_comp_list<-c()
  for (tp_idx in 1:tp_iter){
    
    tp<-unique(timepoints(time_object))[tp_idx]
    if (do_temporal==FALSE){
      return_list<-prep_tp_matrix(time_object,tp,group_names,microarr_dta=rna_biop_dataa)
      comparison_name<-paste0(group_names[2],'_vs_',group_names[1],'_TP_',tp)
    }else{
      tp_next<-unique(timepoints(time_object))[tp_idx+1]
      return_list<-prep_tp_matrix(time_object,tp,c(tp,tp_next),microarr_dta=rna_biop_dataa,do_temporal=TRUE)
      comparison_name<-paste0('TP_',tp_next,'_vs_TP_',tp)
    }
    
    my_mm<-return_list[[1]]
    all_samp_IDs<-return_list[[2]]

    
    my_comp_list<-c(my_comp_list,comparison_name)
    my_eb_res<-calculate_EB(rna_biop_dataa,my_mm,comparison_name)
    
    #Create DESeq2 style results, store total number of genes and significant genes
    temp_gene_num_lst <- produce_DE_results_for_limma(eb_res=my_eb_res,limma_obj = rna_biop_dataa,used_samples = all_samp_IDs,
                                                      filter_choice=filter_choice,save_path = my_save_path,
                                                      gtf_path = gtf_path,gene_int = genes_of_interest,species = species)
    
    #Store number of total and significant genes in list
    gene_num_list[[names(temp_gene_num_lst)]]<-temp_gene_num_lst[[names(temp_gene_num_lst)]]
    
    #Supplement missing gene_ids with systemic names
    my_eb_res$genes$gene_id <- ifelse(my_eb_res$genes$GeneName=='', my_eb_res$genes$SystematicName, my_eb_res$genes$GeneName)
    #Store eb_res results in a list, will be used for fgsea analysis.
    eb_res_list[[comparison_name]]<-my_eb_res
    
  }
  #Save all empirical bayesian results
  save(eb_res_list, file=paste0(my_save_path,'Robjects_checkpoints/final_results.rdata'), version = 2)
  
  #Create a summary table showing total number of DE genes and number of significant DE genes
  create_DE_summary_table(gene_num_list,save_path=my_save_path)
  #Creates the results for the genes of interest
  if (length(genes_of_interest)>0){
    if (is.null(gtf_path)==T){
      target_file<-'DE_raw_data.csv'
    }else{
      target_file<-'DE_results_with_gtf.csv'
    }
    
    create_tables_genes_interest(genes_of_interest,my_comp_list,pval_or_padj=filter_choice,DE_name = target_file,
                                 path_to_DE = my_save_path,
                                 results_path=paste0(my_save_path,'tables_genes_interest/'))
  }
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which retrieves significant DEGs based on given thresholds
#' 
#' The function will search given result files for significant DEGs
#' It requires the path to the differential expression results along
#' with the name of the differential expression results file. It will then search
#' either conditionally or temporally through the file and identify significant
#' genes based on given thresholds. 
#' The results are stored in a list where a merged_df will be all significant genes
#' and their values, the sorted_df will be all unique significant genes sorted by the 
#' filter choice (either adjusted p value or pvalue), the unique_genes will be
#' a list of all significant unique genes.
#' 
#' 
#' @param time_object A time series object
#' @param group_names The conditional names given to the experiment
#' @param DE_path The path to the DE_file
#' @param DE_file_name The name of the file containing the differential expression results.
#' @param p_filter The filter used on 'filter_choice' (either adjusted p value or p value)
#' @param l2fc_filter the log2foldchange filter
#' @param filter_choice The p filter to use (either padj or pvalue)
#' @param do_temporal A boolean indicating if the results are form a temporal comparison
#' or a conditional comparison.
#' 
#' return DEG_list A list containing all significant DEGs, sorted significant DEGs
#' and a vector of unique significant DEGs
#'
#' @examples
#' 
#'
#' @export
find_significant_DEGs <- function(time_object,group_names,DE_path,DE_file_name,p_filter=0.05,l2fc_filter=1,
                                         filter_choice='padj',do_temporal=FALSE){
  
  my_groups<-paste0(group_names[2],'_vs_',group_names[1])
  important_core_columns<-c('gene_id','log2FoldChange','baseMean','pvalue','padj')
  important_col_list<-list()
  
  experiments_to_check<-c()
  if (do_temporal==FALSE){
    for (tp in unique(TimeSeriesExperiment::timepoints(time_object))){
      tp_group<-paste0(my_groups,'_TP_',tp)
      experiments_to_check<-c(experiments_to_check,tp_group)
      samples_used<-timepoints(time_object)
      samples_used<-names(samples_used[samples_used==tp])
      
      important_col_list[[tp_group]]<-c(important_core_columns,samples_used)
    }
  }else{
    my_timepoints=unique(TimeSeriesExperiment::timepoints(time_object))
    for (tp_idx in 1:(length(my_timepoints)-1)){
      tp_group<-paste0('TP_',my_timepoints[tp_idx+1],'_vs_TP_',my_timepoints[tp_idx])
      
      experiments_to_check<-c(experiments_to_check,tp_group)
      
      samples_used_1<-names(which(timepoints(time_object)==my_timepoints[tp_idx]))
      samples_used_2<-names(which(timepoints(time_object)==my_timepoints[tp_idx+1]))
      samples_used<-c(samples_used_1,samples_used_2)
      
      important_col_list[[tp_group]]<-c(important_core_columns,samples_used)
    }
  }

  DEG_list<-list()
  duplicate_gene_names<-c()
  DEG_list[['merged_df']]<-as.data.frame(NULL)
  

  
  for (exp in experiments_to_check){
    print(exp)
    # print(important_col_list[[exp]])
    DE_info<-read_csv(paste0(DE_path,'/',exp,'/',DE_file_name))
    DE_info<-DE_info[DE_info[filter_choice]<p_filter,]
    DE_info<-DE_info[abs(DE_info$log2FoldChange)>=l2fc_filter,]
    
    DE_info_core<-DE_info[,important_core_columns]
    DE_info_samples<-DE_info[important_col_list[[exp]]]
    DEG_list[[exp]]<-DE_info_samples
    
    if (nrow(DEG_list[['merged_df']])==0){
      DEG_list[['merged_df']]<-DE_info_core
    }else{
      print(colnames(DEG_list[['merged_df']]))
      DEG_list[['merged_df']]<-rbind(DEG_list[['merged_df']],DE_info_core)
    }
    duplicate_gene_names<-c(duplicate_gene_names,DE_info_core$gene_id)
  }
  duplicate_gene_names<-duplicate_gene_names[duplicated(duplicate_gene_names)]
  
  #Remove duplicates
  if (filter_choice=='log2FoldChange'){
    DEG_list[['sorted_df']]<-DEG_list[['merged_df']][order(abs(DEG_list[['merged_df']])[[filter_choice]]),]
  }else{
    DEG_list[['sorted_df']]<-DEG_list[['merged_df']][order(DEG_list[['merged_df']][[filter_choice]]),]
  }
  DEG_list[['sorted_df']]<-DEG_list[['sorted_df']][!duplicated(DEG_list[['sorted_df']]$gene_id), ]
  
  DEG_list[['unique_genes']]<-DEG_list[['sorted_df']]$gene_id
  return(DEG_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates a combined list of both temporal and conditional DEGS
#' 
#' The function takes in two DEG_lists as created by \code{find_significant_DEGs}
#' one representing significant temporal DEGs and the other the conditional 
#' significant DEGs.
#' The function will create a dataframe with the unique genes between the 
#' two lists
#' 
#' 
#' @param deg_list_temporal The DEG list for temporally significant genes
#' @param deg_list_conditional The DEG list for the conditionally significant genes
#'
#'return merged_df A dataframe containing the merged results
#'
#' @examples
#' 
#'
#' @export
merge_conditional_temporal_df<-function(deg_list_temporal,deg_list_conditional){
  conditional_df<-as.data.frame(deg_list_conditional['sorted_df']$sorted_df)
  conditional_df<-conditional_df[,c('gene_id','padj')]
  colnames(conditional_df)=c('gene_id','padj_cond')
  temporal_df<-as.data.frame(deg_list_temporal['sorted_df']$sorted_df)
  temporal_df<-temporal_df[,c('gene_id','padj')]
  colnames(temporal_df)=c('gene_id','padj_temp')
  
  merged_df<-merge(conditional_df,temporal_df,all=TRUE)
  min_vect<-c()
  for (idx in 1:nrow(merged_df)){
    min_vect<-c(min_vect,min(merged_df[idx,2:3],na.rm = TRUE))
  }
  merged_df$min_padj<-min_vect
  merged_df<-merged_df[order(merged_df[['min_padj']]),]
  
  return(merged_df)
}



# Custom heatmap functions -------------------

#' A wrapper function which calls the necessary functions to create a heatmap
#' illustrating the significant DEGs of a given list originating from 
#' \code{find_significant_DEGs}
#' The function will automatically save the file to the relevant location (TS_results).
#' 
#' @param time_object A time object
#' @param group_names The groups used in the time series object
#' @param DEG_list a list of significant degs as created by \code{find_significant_DEGs}
#' @param log_transform Boolean indicating if the results should be log transformed
#' on the heatmap
#' @param plot_file_name The file name to be given to the heatmap once it is saved

#'
#' @examples
#' 
#'
#' @export
custom_heatmap_wrapper<-function(time_object,group_names,DEG_list,log_transform=T,plot_file_name='custom_DEG_heatmap'){
  my_heat_list<-create_heatmap_matrix(time_object,group_names=group_names,DEG_list)

  
  temp_list<-prepare_heat_data(my_heat_list,log_transform = log_transform)
  
  my_heat_mat<-temp_list[['heat_matrix']]
  my_region_split<-temp_list[['region_split']]
  my_group_split<-temp_list[['group_split']]
  my_l2fc_vect<-temp_list[['l2fc_vector']]
  
  if (log_transform==T){
    my_l2fc_vect<-log_transform_l2fc_vect(my_l2fc_vect)
  }
  
  
  plot_custom_DE_heatmap(my_heat_mat,my_region_split,my_group_split,my_l2fc_vect,log_transform = log_transform,
                         save_path='TS_results/',plot_file_name = plot_file_name)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates a data matrix to be used in the creation of a custom heamap
#' 
#' This function is intended to be used within the \code{cusom_heatmap_wrapper}
#' 
#' The function creates a list containing the main data matrix along with three vectors
#' l2fc_vector: A named vector containing the L2FC value for each gene of the matrix
#' group_vector: A vector containing the new names of the groups to represent the group
#' they are in
#' gene_vector: A named vector containing the gene names along with the group that
#' they are in.
#' 
#' The vectors are created to account for the possibility of duplicates. For example, 
#' a gene may be significant in two timepoint comparisons, and therefore have two
#' different values. The vectors give the genes unique names based on their grouping
#' in order to enable them to have different values if necessary.
#' 
#' 
#' @param time_object A time series object
#' @param group_names The names given to the conditionals of the time series object
#' @param DEG_list a list as produced by \code{find_significant_DEGs}
#' 
#' return heat_list A list contianing the main data matrix (count values) and the 
#' three vectors described above.
#' 
#' 
#' @examples
#' 
#'
#' @export
create_heatmap_matrix<-function(time_object,group_names,DEG_list){
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
  
  #Iterate over experiments/timepoints
  for (tp_group in exps_interest){
    # tp_group<-paste0(my_groups,'_TP_',tp)#Create name to look for
    temp_df<-as.data.frame(DEG_list[[tp_group]])#Extract df associated to timepoint

    genes_mod<-temp_df$gene_id#Extract gene names
    genes_mod<-paste0(genes_mod,'_TP',tp_group)#Modify gene names by adding time point

    #Extract and name log2FoldChange
    temp_l2fc<-temp_df$log2FoldChange
    names(temp_l2fc)<-genes_mod
    
    temp_df<-temp_df[,6:ncol(temp_df)]#Select for experiment related data

    row.names(temp_df)=genes_mod
    temp_df<-t(temp_df)#transpose
    row.names(temp_df)=unname(replicates(time_object)[row.names(temp_df)])#rename using replicates
    

    
    #Fill the main dataframe
    if (nrow(heat_df)==0){
      heat_df<-temp_df
    }else{
      heat_df<-cbind(heat_df,temp_df)
    }
    #Concatenate the log2FoldChanges
    l2fc_vect<-c(l2fc_vect,temp_l2fc)
    temp_gene_vector<-rep(paste0(tp_group,' (',length(genes_mod),')'),length(genes_mod))
    gene_vector<-c(gene_vector,temp_gene_vector)
  }
  
  my_group_df<-as.data.frame(TimeSeriesExperiment::groups(time_object))
  my_group_df<-cbind(my_group_df,unname(replicates(time_object)))
  colnames(my_group_df)<-c('Group','Replicate')
  my_group_df<-unique(my_group_df)
  row.names(my_group_df)=my_group_df$Replicate
  
  group_vector<-my_group_df[row.names(heat_df),]
  group_vector<-group_vector$Group

  #Fill and return result list
  heat_list[['main_matrix']]<-heat_df
  heat_list[['l2fc_vector']]<-l2fc_vect
  heat_list[['group_vector']]<-group_vector
  heat_list[['gene_vector']]<-gene_vector
  return(heat_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which adjusts and renames certain elements of the heat_list
#' produced by \code{create_heatmap_matrix}
#' 
#' It performs a log_transformation if necessary and converts the vectors into 
#' factors to be used for annotation purposes within the heatmap
#' 
#' This function is intended to be used within the \code{cusom_heatmap_wrapper}
#' 
#' @param matrix_list The heat_list produced by \code{create_heatmap_matrix}
#' @param log_transform If a log transformation of the data is necessary
#'
#' return return_list An updated version of the inputted matrix_list
#'
#' @examples
#' 
#'
#' @export
prepare_heat_data <- function(matrix_list,log_transform){
  heat_matrix<-matrix_list[['main_matrix']]
  
  if (log_transform==T){
    heat_matrix<-log10(as.matrix(heat_matrix+1))
  }
  
  
  group_split<-matrix_list[['group_vector']]

  numbered_regions<-matrix_list[['gene_vector']]

  #Sets up the order of appearance
  region_split <- factor(numbered_regions, levels=unique(numbered_regions))
  group_split <- factor(group_split, levels=unique(group_split))

  l2fc_vector<-matrix_list[['l2fc_vector']]
  names(l2fc_vector)=region_split
  return_list<-list('heat_matrix'=heat_matrix,'region_split'=region_split,
                    'group_split'=group_split,'l2fc_vector'=l2fc_vector)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which log transforms a vector containing log2FoldChange values
#' 
#' The function first identifies the location of negative values. It then performs
#' the log transformation on the absolute values of the vector.
#' Lastly it returns originally negative values to negative values
#' This is done as log transformation does not work on negative values, therefore a
#' work around had to be done.
#' 
#' This function is intended to be used within \code{custom_heatmap_wrapper}
#' 
#' @param l2fc_vector The log2FoldChange vector to be log transformed
#'
#' return new_l2fc The log transformed vector
#'
#' @examples
#' 
#'
#' @export

log_transform_l2fc_vect <-function(l2fc_vector){
  neg_idx<-c()
  for (idx in 1:length(l2fc_vector)){
    if (l2fc_vector[idx]<0){
      neg_idx<-c(neg_idx,idx)
    }
  }
  
  
  new_l2fc<-log10(abs(l2fc_vector))
  for (idx in 1:length(new_l2fc)){
    if (idx %in% neg_idx){
      new_l2fc[idx]<- -new_l2fc[idx]
    }
  }
  
  return(new_l2fc)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which plots a custom heatmap using the complexHeatmap package
#' 
#' The function is intended to be called within the \code{custom_heatmap_wrapper}
#' 
#' The function creates a segmented heatmap where rows are samples/patients and columns
#' are genes. The heatmap is column segmented based on the number of differential expression
#' analyses contained within it. A log2foldChange histogram is attached to the bottom
#' of the heatmap to indicate significance of the genes being visualized.
#' 
#' The column/segment legend is given at the bottom of the heatmap with the number
#' of genes in each group being also provided.
#' 
#' The gradient legend of the values within the heatmap is given on the right hand side.
#' 
#' The heatmap is saved in both png and svg format in the declared location via the
#' 'save_path' parameter.
#' 
#' 
#' @param heat_mat The main count matrix of the heatmap
#' @param col_split The named vector used to split the columns into it's segments
#' @param row_splits The named vector used to split the rows into the two groups
#' @param l2fc_col The log2foldchange vector used to create the histogram
#' @param log_transform boolean indicating if the results were log transformed or not
#' @param save_path The path to which the heatmap will be saved
#' @param plot_file_name The name given to the saved heatmap
#' @param custom_width The width of the heatmap
#' @param custom_height The height of the heatmap
#'
#' @examples
#' 
#'
#' @export
plot_custom_DE_heatmap <-function(heat_mat,col_split,row_splits,l2fc_col, log_transform,
                                  save_path='',plot_file_name='custom_heatmap',
                                  custom_width=15,custom_height=5){
  
  if(log_transform==T){
    bottom_histo = HeatmapAnnotation('log10(l2FC)' = anno_barplot(l2fc_col))
    count_legend<-'log10(counts)'
  }else{
    bottom_histo = HeatmapAnnotation('l2FC' = anno_barplot(l2fc_col))
    count_legend<-'counts'
  }
  
  fill_set_regions=gpar(fill=2:(length(unique(col_split))+1))
  top_annot = HeatmapAnnotation(foo = anno_block(gp = fill_set_regions, labels = rep(NULL,length(unique(col_split)))))
  
  fill_set_groups=gpar(fill = c("#998ec3", "#f1a340"))
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
                cluster_column_slices=F,cluster_row_slices = F),
       annotation_legend_list = lgd,
       annotation_legend_side = 'bottom'
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
                cluster_column_slices=F,cluster_row_slices = F),
        annotation_legend_list = lgd,
        annotation_legend_side = 'bottom'
  )
  dev.off()
}

# PART processing and heatmap  -------------------

#' Function which formats the count values for a PART heatmap
#' 
#' This function is intended to be used within \code{PART_heat_map}
#' 
#' The function retrieves the counts from the time series object and filters
#' it to only contain relevant genes and samples. The data can also be scaled
#' if the parameter 'scale' is set to T
#' 
#' 
#' @param object A time series object
#' @param target_genes Vector of genes to keep/use
#' @param scale Boolean indicating if the data should be scaled
#' @param use_samples Vector of samples to keep/use
#' 
#' return Y The count matrix for relevant genes and samples
#'
#' @examples
#' 
#'
#' @export
prep_counts_for_PART_heat <-function(object,target_genes,scale,use_samples){
  #Retrieve and prep assay
  cnts <- assays(object)$norm
  rownames(cnts) <- rowData(object)[, 'symbol']
  colnames(cnts) <- colData(object)[, 'sample']
  
  found_groups<-setNames(object@group,colnames(object))
  found_replicates<-setNames(object@replicate,colnames(object))
  found_timepoints<-setNames(object@timepoint,colnames(object))
  
  if (is.null(use_samples)==F){
    cnts <- cnts[, use_samples]
    found_groups =found_groups[use_samples]
    found_replicates =found_replicates[use_samples]
    found_timepoints =found_timepoints[use_samples]
  }
  
  if (is.null(target_genes)==T){
    top_feat <- apply(cnts, 1, sd)
    top_feat <- names(top_feat)[order(-top_feat)[seq_len(num.feat)]]
  }else{
    top_feat <- target_genes
    gene_not_found <- top_feat[!top_feat %in% rownames(cnts)]
  }
  
  
  top_feat <- top_feat[top_feat %in% rownames(cnts)]
  
  cols_ordered <- order(found_groups, found_replicates, 
                        found_timepoints)
  
  #creates, orders, and scales the normalized count matrix for the top features selected
  Y <- cnts[top_feat, cols_ordered]
  if (scale) {
    Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
  }
  
  return (Y)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which computes the PART clusters for a given matrix
#' 
#' This function is intended to be used within \code{PART_heat_map}
#' 
#' The function set's a see for reproducibility of results, it then uses the
#' part function from the clusterGenomics package to establish which genes belong
#' to what clusters.
#' The function then retrieves this information and creates the rowAnnotation
#' object that will be used in the complexHeatmap to illustrate the results
#' 
#' The function may also save the PART results as a csv if indicated to do so
#' 
#' 
#' @param main_matrix The matrix used for PART clustering
#' @param save_name the file name to save the PART data, if NULL no saving
#' @param part_recursion The number of recursions for PART calculation
#' @param part_min_clust The minimum number of genes per cluster
#' @param dist_param The distance parameter for clustering
#' @param hclust_param The hierarchical clustering method/parameter to be used 
#'
#' return list(row_annot,rowclust) row_annot is the rowAnnotation object for the 
#' complexHeatmap, rowClust is the hierarchical clustering results.
#'
#' @examples
#' 
#'
#' @export
compute_PART<-function(main_matrix,save_name='PART_data.csv',part_recursion=50,part_min_clust=10,
                       dist_param="euclidean", hclust_param="average"){
  #Calculates the clustering using the 'part' algorithm
  set.seed('123456')
  calculated_clusters = part(main_matrix,B=part_recursion,minSize=part_min_clust)
  rowclust = hclust(dist(main_matrix,method=dist_param),method=hclust_param)
  
  
  #Sets up the format for the color bar which will illustrate the different clusters found
  num_clusts <- length(unique(calculated_clusters$lab.hatK))
  cols <- RColorBrewer::brewer.pal(8, name = "Set3")[seq_len(min(8, num_clusts))]
  clust_cols <- grDevices::colorRampPalette(colors = cols)(num_clusts)
  names(clust_cols) <- unique(as.character(calculated_clusters$lab.hatK))
  
  clust_ordered <- unique(as.character(calculated_clusters$lab.hatK[rowclust$order]))
  
  #Rename clusters so they appear in order (i.e. 1,2,3,4)
  my_iter<-1
  new_vect<-c()
  for(val in unique(calculated_clusters$lab.hatK[rowclust$order])){
    rep_value<-unname(table(calculated_clusters$lab.hatK[rowclust$order])[as.character(val)])
    new_vect<-c(new_vect,rep(as.character(my_iter),rep_value))
    my_iter<-my_iter+1
  }
  #Cluster illustration stored as rowAnnotation
  row_annot <- ComplexHeatmap::rowAnnotation(gene_cluster = as.character(calculated_clusters$lab.hatK), 
                                             col = list(gene_cluster=clust_cols),
                                             show_annotation_name=F,
                                             annotation_legend_param = list(title = "clusters", at = clust_ordered,
                                                                            labels = unique(new_vect)))
  
  
  #Create the color vector to add to file
  col_vect<-clust_cols[clust_ordered]
  names(col_vect)=seq(1,length(col_vect))
  color_vector<-c()

  #Saves the file as a csv if save_name is not null
  if (is.null(save_name)==F){
    FD_2 <- main_matrix[rowclust$order,]
    FD_2 <- as.data.frame(FD_2)
    FD_2 <- add_column(FD_2, gene_cluster = new_vect, .before = 1)
    for (clust in unique(FD_2$gene_cluster)){
      color_vector<-c(color_vector,rep(col_vect[clust],nrow(FD_2[FD_2$gene_cluster==clust,])))
    }
    FD_2 <- add_column(FD_2, cluster_col = color_vector, .after = 1)
    write.csv(FD_2,save_name, row.names = T)
  }
  
  return_list<-list(row_annot,rowclust)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which prepares the top annotation for the PART heatmap
#' 
#' This function is intended to be used within \code{PART_heat_map}
#' 
#' The top annotation shows annotation blocks for the groups as well as annotation
#' gradients for the time points.
#' Two sets of top annotations are created, one labelled and one unlabeled.
#' In addition, the column splits are created based on the replicates
#' 
#' 
#' @param object A time series object
#' @param main_matrix The value matrix that will be used in the PART heatmap
#' 
#' return list(top_annot_labels, top_annot_no_labels,col_split)
#' The two types of top annotaiton and the named vector used to create the column splits
#'
#' @examples
#' 
#'
#' @export
prepare_top_annotation_PART_heat<-function(object,main_matrix){
  found_timepoints<-TimeSeriesExperiment::timepoints(object)[colnames(main_matrix)]
  found_replicates<-TimeSeriesExperiment::replicates(object)[colnames(main_matrix)]
  
  col_split<-factor(unname(found_replicates),levels=unique(unname(found_replicates)))
  
  num_cols<-table(sapply(strsplit(levels(col_split),"_"), `[`, 1))
  fill_set_groups=gpar(fill=c(rep("#1f78b4",unname(num_cols[1])),rep("#e31a1c",unname(num_cols[2]))))
  
  my_cols<-c('#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
  # If possible, use custom colors, otherwise use colorbrewer
  if (length(unique(found_timepoints))>length(my_cols)){
    num_tp <- length(unique(found_timepoints))
    cols <- RColorBrewer::brewer.pal(8, name = "Set1")[seq_len(min(8, num_tp))]
    time_cols <- grDevices::colorRampPalette(colors = cols)(num_tp)
    names(time_cols) <- unique(found_timepoints)
  }else{
    time_cols<-my_cols[1:length(unique(found_timepoints))]
    names(time_cols)<-unique(found_timepoints)
  }
  
  
  top_annot_labels = HeatmapAnnotation(foo = anno_block(gp = fill_set_groups, labels = unique(col_split)),
                                       timepoints=found_timepoints,
                                       col = list(timepoints = time_cols),
                                       show_annotation_name=F,
                                       annotation_legend_param = list(title = "timepoints", at = unique(found_timepoints),
                                                                      labels = unique(found_timepoints)))
  
  top_annot_no_labels = HeatmapAnnotation(foo = anno_block(gp = fill_set_groups, labels = rep(NULL,length(unique(col_split)))),
                                          timepoints=found_timepoints,
                                          col = list(timepoints = time_cols),
                                          show_annotation_name=F,
                                          annotation_legend_param = list(title = "timepoints", at = unique(found_timepoints),
                                                                         labels = unique(found_timepoints)))
  
  return(list(top_annot_labels,top_annot_no_labels,col_split))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A wrapper function used to create a heatmap illustrating the results from PART
#' clustering
#' 
#' The heatmap shows selected genes as rows and replicates of the time series experiment
#' as columns.
#' Gene clusters are identified by colored row annotation. Replicates are split
#' by their grouping (conditional) and ordered (within each other) via time points
#' The legend for groupings, clusters, timepoints, and heatmap values is given on
#' the right hand side of the heatmap.
#' 
#' The heatmap is saved in the TS_results folder. The heatmap is saved twice, once 
#' in png format, and the other in svg format
#' 
#' 
#' @param object A time series object
#' @param target_genes The genes to cluster and plot
#' @param group_graph_vector A named vector showing the groups used along with the desired color
#' to label them with
#' @param heat_name The file name given to the saved heatmap
#' @param scale Boolean indicating if the data should be scaled
#' @param use_samples Indication of which samples to use, can be NULL (all samples used)
#'
#' @examples
#' 
#'
#' @export
PART_heat_map<-function(object, target_genes, group_graph_vector, heat_name='custom_heat_map',scale=T, use_samples=NULL){
  
  
  Y<-prep_counts_for_PART_heat(object,target_genes,scale,use_samples)
  
  PART_save<-paste0(heat_name,'_data.csv')
  part_results<-compute_PART(Y,save_name = PART_save,part_min_clust = 30)
  row_annot<-part_results[[1]]
  rowclust<-part_results[[2]]
  
  
  top_annot_results<-prepare_top_annotation_PART_heat(object,Y)
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
  lgd = Legend(labels = names(group_graph_vector), title = "groups", legend_gp = gpar(fill = unname(group_graph_vector)))
  
  
  #Plot the heatmap
  svg(paste0(heat_name,"_with_names.svg"),height=30,width=20)
  draw(
    ComplexHeatmap::Heatmap(
      Y, name = "Z-score", cluster_columns = F,cluster_rows=rowclust, show_column_dend = T,
      row_names_gp = grid::gpar(fontsize = 8),left_annotation = row_annot,
      show_row_names = T,top_annotation = top_annot_labels,column_split = col_split,cluster_column_slices = T,
      column_gap = unit(gap_vect, "mm"),show_column_names = T,border=F,column_title = NULL),
    annotation_legend_list = lgd
  )
  dev.off()
  
  
  svg(paste0(heat_name,".svg"))
  draw(
    ComplexHeatmap::Heatmap(
      Y, name = "Z-score", cluster_columns = F,cluster_rows=rowclust, show_column_dend = T,
      row_names_gp = grid::gpar(fontsize = 8),left_annotation = row_annot,
      show_row_names = F,top_annotation = top_annot_no_labels,column_split = col_split,cluster_column_slices = T,
      column_gap = unit(gap_vect, "mm"),show_column_names = F,border=F,column_title = NULL),
    annotation_legend_list = lgd
  )
  dev.off()
}

# plotting functions -------------------


#' A wrapper function which plots the cluster trajectory of all clusters within
#' a cluster map. 
#' 
#' The function limits the number of cluster per plot to 20 clusters
#' as to not overcrowd the figure. Supplementary figures will be created to ensure
#' that all clusters will be plotted. 
#' Cluster trajectories are saved as svg files
#' 
#' @param time_object A time series object
#' @param group_names The group names given to the time series object and shown
#' on the plots
#' @param clust_map The cluster map
#' @param graphic_vector A named vector stating group names and their associated color
#' @param plot_name The location and name to which the plots/figures will be saved
#' @param log_timepoints Boolean indicating if the timepoints should be logged. Can 
#' be practical if timepoints are not at even intervals.
#'
#' @examples
#' 
#'
#' @export
my_clust_traj_wrapper<-function(time_object,group_names,clust_map,graphic_vector,
                                plot_name='TS_results/cluster_traj',log_timepoints=F){
  
  clust_map<-clust_map[,c(1,2)]
  #Create and sort table for number of clusters and genes per cluster
  clust_overview<-table(clust_map$cluster)
  clust_overview<-clust_overview[order(-clust_overview)]
  
  #10 cluster per figure is optimal
  num_needed_figures<-ceiling(length(clust_overview)/10)
  
  #Iterate over number of necessary figures
  for (idx in 1:num_needed_figures){
    #Filter the cluster_map dataframe for the required clusters
    max_clust=10*idx
    if (idx==1){
      min_clust<-1
    }else{
      min_clust<-10*(idx-1)+1
    }
    clusters_to_plot<-names(clust_overview)[min_clust:max_clust]
    clusters_to_plot<-clusters_to_plot[!is.na(clusters_to_plot)]#Remove NAs
    
    sub_clust_map<-clust_map[clust_map$cluster %in% clusters_to_plot,]
    
    if (num_needed_figures > 1){
      save_name<-paste0(plot_name,'_',idx,'_of_',num_needed_figures,'.svg')
    }else{
      save_name<-paste0(plot_name,'.svg')
    }
    
    svg(save_name,width=12,height=12)
    print(plotTS_clust_mod(c_map=sub_clust_map,object=time_object, groups.labs=group_names,transparency=0.2,scales = 'free_x',log_timepoints = log_timepoints)+
            scale_color_manual(values = graphic_vector) +
            theme(strip.text = element_text(size = 10)))
    dev.off()
  }
}







### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





#' A function which will plot the time series based on a custom cluster_map
#' 
#' 
#' The function takes in a cluster map, which is a file containing gene names and
#' their associated cluster. This can be obtained by the \code{create_cluster_map_from_csv}
#' function.
#' 
#' For a single gene, the function will create plots showing the genes expression
#' over the time points for each group. The solid line of the plot shows the mean
#' expression of all the replicates put together, the lighter areas show the actual
#' expression per replicate.
#' 
#' This plot effectively allows the visualization of the expression of genes from a 
#' specific cluster of the cluster_map.
#' It is practical to use in conjunction with a heatmap to further explore a cluster 
#' found in the heatmap.
#' 
#' 
#' @param c_map The cluster map
#' @param object The time series object
#' @param groups.labs The different groups contained in the time series object to be plotted
#' @param transparency transparency of trajectory lines.
#' @param ncol number of columns in the factet plot.
#' @param scales character scalar indecating facet scales, by default "free".
#'
#' @return plt The plot itself
#' 
#' @examples
#' plotHeatmap_part_mod(TS_object, num.feat = 300, save_name = 'TS_results/top_300_variable_genes.csv')
#' 
#' cluster_map <- create_cluster_map_from_csv(TS_object,csv_path='TS_results/top_300_variable_genes.csv')
#' 
#' plotTS_clust_mod(c_map=cluster_map,object=TS_object, groups.labs=c('WT','Loxp'),transparency=0.2)+
#'   scale_color_manual(values = c("WT" = "#1f78b4", "Loxp" = "#e31a1c")) +
#'   theme(strip.text = element_text(size = 10)) +
#'   ylim(NA, 0.55)
#' @export
plotTS_clust_mod <- function (c_map, object, groups.labs,transparency = 0.5, ncol = 4, scales = "free",
                              log_timepoints=F){
  
  cluster_map <-c_map
  
  #freq_df contains the number of genes in each cluster, created for visualiztion on the plots
  freq_df <- cluster_map %>%
    group_by(cluster) %>%
    summarise(freq = n()) %>%
    arrange(desc(freq))
  
  #Using the freq df, creates the format for the illustration of the amount of
  #genes in each cluster
  cat_levels <- paste0("[", freq_df$cluster, ": ", freq_df$freq, "]")
  cat_levels <- paste(rep(cat_levels, each = length(groups.labs)),
                      rep(groups.labs, length(cat_levels)))

  #Stores the necessary features and creates ts_data, which contains all the
  #information to be plotted in the trajectory plots
  features <- cluster_map$feature
  
  # cluster_map$feature <- cluster_map$symbol
  # cluster_map<-cluster_map[c('feature','cluster')]
  
  ts_data <- 
    timeSeries(object, "ts_collapsed") %>%
      filter(feature %in% features) %>%
      dplyr::select(-replicate, -contains("Lag_")) %>%
      left_join(cluster_map) %>%
      gather(
        key = "timepoint", value  = "value",
        -feature, -group, -cluster
      ) %>%
      left_join(freq_df) %>%
      mutate(
        timepoint = as.numeric(timepoint),
        value = as.numeric(value),
        category = paste0("[", cluster, ": ", freq, "] ", group)
      ) %>%
      arrange(cluster, group, feature, timepoint) %>%
      mutate(
        category = factor(category, levels = cat_levels))
  View(ts_data)
  # Compute cluster mean expression profile for each each gene in each group
  ts_cluster_mean <-
    ts_data %>%
      dplyr::select(-feature) %>%
      group_by(cluster, group, category, timepoint) %>%
      summarize_all(mean) %>%
      left_join(freq_df) %>%
      arrange(cluster, group)
  ts_data$group <- as.character(ts_data$group)

  #Creates the actual plot
  min_y_log <- min(ts_data$value)
  max_y_log <- max(ts_data$value)
  if (log_timepoints==T){
    plt <- ggplot(ts_data, aes(y = value , x = log10(timepoint), color = group))+
      geom_line(aes(group = feature), alpha = transparency) +
      geom_point() +
      geom_line(
        data = ts_cluster_mean, lwd = 1.5, color = "grey50",
        aes(group = group)
      ) +
      # ylim(min_y_log,max_y_log) +
      facet_wrap(~category, scales = scales, ncol = ncol)
  }else{
    plt <- ggplot(ts_data, aes(y = value , x = timepoint, color = group))+
      geom_line(aes(group = feature), alpha = transparency) +
      geom_point() +
      geom_line(
        data = ts_cluster_mean, lwd = 1.5, color = "grey50",
        aes(group = group)
      ) +
      # ylim(min_y_log,max_y_log) +
      facet_wrap(~category, scales = scales, ncol = ncol)
  }

  
  return(plt)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' An adaption of the timeseries plotting function to allow us to plot
#' specific genes from clusters.
#' 
#' 
#' @param object timeSeries objects
#' @param features the features to be plotted
#' @param timepoints_selected The timepoints to be plotted
#' @param trans
#' @param smooth Weather or not the trajectory should be smoothed
#' @param ncol
#' @param scales the type of scaling to use
#' @param smooth_shadow If shadowing should be used on the smoothing
#'
#' @examples
#' 
#'
#' @export
plotTimeSeries_mod <- function(object, features = rownames(object), 
                               timepoints_selected=NULL,
                           trans = FALSE, smooth = TRUE, ncol = 5,
                           scales = "free",smooth_shadow=T)
{
  feature <- symbol <- timepoint <- value <- group <- category <- NULL
  if (!is(object, "TimeSeriesExperiment")) 
    stop("Input must be a 'TimeSeriesExperiment' object.")
  if (!validObject(object))
    stop("Invalid TimeSeriesExperiment object.")
  if(!all(features %in% rownames(object)))
    stop("'features' must be a subset of rownames(object)")
  
  if (!"ts" %in% names(timeSeries(object))) {
    object <- makeTimeSeries(object)
  }
  feature_data <- rowData(object) %>% 
    as.data.frame() %>%
    filter(feature %in% features) %>%
    arrange(factor(feature, levels = features))
  
  if(!"symbol" %in% colnames(feature_data)){
    feature_data$symbol <- feature_data$feature
  }
  
  if(trans) {
    ts_data <- timeSeries(object, "ts_trans")
  } else {
    ts_data <- timeSeries(object, "ts")
  }
  ts_data <- suppressMessages(
    ts_data %>%
      filter(feature %in% features) %>%
      dplyr::select(-starts_with("Lag_")) %>%
      gather(key = "timepoint", value = "value", -(feature:replicate)) %>%
      left_join(feature_data %>% dplyr::select(feature, symbol)) %>%
      mutate(
        symbol = factor(symbol, levels = feature_data$symbol),
        timepoint = as.numeric(timepoint),
        category = paste0(group, "_", replicate)) 
  )
  
  if (is.null(timepoints_selected)==F){
    ts_data <- ts_data[ts_data$timepoint %in% timepoints_selected,]
  }else{
    found_zeros<-which(ts_data$timepoint==0)
    ts_data$timepoint<-ts_data$timepoint
    ts_data$timepoint<-	replace(ts_data$timepoint, found_zeros, 0)
  }
  
  # View(ts_data)
  
  plt <- ggplot(
    ts_data,
    aes(x = timepoint, y = value, color = group)) +
    geom_point(size = 1) +
    facet_wrap(~ feature, scales = scales, ncol = ncol)
  
  if(length(unique(ts_data$replicate)) > 1) {
    plt <- plt + geom_line(aes(group = category), lty = 3, alpha = 0.7)
  }
  if(smooth) {
    plt <- plt + geom_smooth(aes(x = timepoint ), lwd = 1.5,se=smooth_shadow)
  } else {
    ts_data_mean <- suppressMessages(
      ts_data %>%
        dplyr::select(-replicate) %>%
        group_by(feature, group, timepoint ) %>%
        summarise(value = mean(value)) %>%
        left_join(feature_data)
    )
    plt <- plt + geom_line(data = ts_data_mean, lwd = 1.5)
  }
  return(plt)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Wrapper function which plots various PCA types
#' 
#' Function will eventually be replaced by a single PCA showing all the information
#' illustrated by the three different PCAs
#' 
#' 
#' @param object A time series object
#' @param graph_vect A named vector stating group names and their associated color
#' 
#' @examples
#' 
#'
#' @export
wrapper_plot_PCAs_TS <- function(object,graph_vect){
  dir.create('TS_results/PCAs')
  #Generates PC components
  TS_object <- runPCA(TS_object, var.stabilize.method = "asinh")
  #Plot PCA for groups
  svg("TS_results/PCAs/group_PCA.svg")
  plotSamplePCA(TS_object, col.var = "group", size = 5)
  dev.off()
  #Plot PCA for timepoints
  svg("TS_results/PCAs/timepoint_PCA.svg")
  plotSamplePCA(TS_object, col.var = "timepoint", size = 5)
  dev.off()
  #PCA showing gene profiles
  svg("TS_results/PCAs/TS_PCA.svg")
  plotTimeSeriesPCA(
    TS_object,
    linecol = graph_vect,
    main = paste("Visualizing gene profiles with PCA"),
    m  = 15, n = 15, col = adjustcolor("grey77", alpha=0.7),
    cex.main = 3, cex.axis = 2, cex.lab = 2)
  dev.off()
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Wrapper function for plotting of gene trajectories
#' 
#' The function plots the gene trajectories of all the genes within a designated cluster
#' It limits the number of genes per figure to 10. It will create the necessary amount
#' of figures to plot all the genes of the cluster
#' 
#' The figures are saved in svg format
#' 
#' @param time_object A time series object
#' @param cluster_map The cluster map
#' @param target_cluster The name of a cluster within the cluster map
#' @param save_location The location where the plots will be saved
#'
#' @examples
#' 
#'
#' @export
wrapper_function_gene_traj <- function(time_object,gene_vector,save_location='TS_results',graph_vect,selected_timepoints=NULL){
  #Keep only genes present in time_object
  gene_vector<-gene_vector[gene_vector %in% rowData(time_object)$symbol]

  num_figs<-ceiling(length(gene_vector)/10)
  
  for (idx in 1:num_figs){
    if (idx == 1){
      min_gene<-1
    }else{
      min_gene<-10*(idx-1)+1
    }
    max_gene<-10*idx
    genes_to_plot<-gene_vector[min_gene:max_gene]
    genes_to_plot<-genes_to_plot[!is.na(genes_to_plot)]#Remove NAs
    
    save_name<-paste0(save_location,'/',idx,'_of_',num_figs,'.svg')
    
    store_res<-plotTimeSeries_mod(time_object, features = genes_to_plot,smooth_shadow = F,scales = 'free_x',timepoints_selected = selected_timepoints) +
      scale_color_manual(values = graph_vect)
    svg(save_name,width=13)
    print(store_res)
    dev.off()
  }
}


# Unused/old functions  -------------------

#' Function which obtains DEGs for a specific time point
#' 
#' The function takes in a DEG_list as created by \code{find_significant_DEGs}
#' The function retrieves the DEGs for the time points of interest as well
#' as the inputted threshold for the adjusted p value and the log 2 fold change
#' 
#' 
#' @param DE_list A list containing the diffferential expression results for each time point
#' @param times_of_interest Vector containing the time points of interest
#' @param adjPval_thresh The adjusted p value threshold to use to identify DEGs
#' @param L2FC_thresh The log2FoldChange threshold to use to identify DEGs
#'
#' 
#' @return DEG_return_list A list splitting the DEGs in three groups- all, upregulated, downregulated
#'
#' @examples
#' plotHeatmap_part_mod(TS_object, num.feat = 300, save_name = 'TS_results/top_300_variable_genes.csv')
#' cluster_map <- create_cluster_map_from_csv(TS_object,csv_path='TS_results/top_300_variable_genes.csv')
#'
#' @export
extract_DEGs_from_time <- function(DE_list,times_of_interest,adjPval_thresh=0.05,L2FC_thresh=1){
  
  for (time in times_of_interest){
    time <- as.character(time)
    DE_extract <- DE_list[[time]][DE_list[[time]]$adj.P.Val < adjPval_thresh,]
    time_up <- DE_extract[DE_extract$logFC > L2FC_thresh,]
    time_down <- DE_extract[DE_extract$logFC < (-L2FC_thresh),]
    
    if (exists('DEG_up')==T){
      DEG_up <- intersect(DEG_up,time_up$symbol)
    }else{
      DEG_up <- time_up$symbol
    }
    
    if (exists('DEG_down')==T){
      DEG_down <- intersect(DEG_down,time_down$symbol)
    }else{
      DEG_down <- time_down$symbol
    }
    
  }
  DEG_return_list <- list()
  DEG_return_list[['upregulated']]<-DEG_up
  DEG_return_list[['downregulated']]<-DEG_down
  DEG_return_list[['all']]<-c(DEG_up,DEG_down)
  
  return(DEG_return_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which extracts the genes involved in a or several GO pathways
#' 
#' The function requires annotation dbi object and string.
#' The annotation object is used to create a godata dataframe. This dataframe is then
#' queried for the GO's submitted in the GO_vector, and converted from entrezID to symbol
#' using the mapIDs function along with the annotaiton dbi string.
#' 
#' 
#' @param organism_obj A annotation dbi object
#' @param organism_str A string version of the annotation dbi object
#' @param GO_vector A vector containing the GO ID(s) of interest
#' @param ontology The ontology to be used
#'
#' 
#' @return GO_gene_list A list which splits the found gene symbols by the GO term used
#' to find them
#'
#' @examples
#' 
#' @export
find_genes_in_GOs <- function(organism_obj, organism_str, GO_vector, ontology='BP'){
  d <- godata(organism_str, ont=ontology)
  #Check if all GOs targeted are in database
  GO_gene_list <- list()
  for (single_GO in GO_vector){
    sub_db <- d@geneAnno[d@geneAnno$GO==single_GO,]
    GO_gene_list[[single_GO]] <- as.vector(mapIds(organism_obj, keys = sub_db$ENTREZID, keytype = "ENTREZID", column="SYMBOL"))
    # GO_gene_list[[single_GO]] <- as.vector(sub_db$EntrezID)
  }
  return(GO_gene_list)
}



clusterTimeSeries_mod <- function(object, n.top.feat = 1000,
                                  groups.selected = "all", lambda = c(0.5, 0.25),
                                  clust.params = list())
{
  if (!is(object, "TimeSeriesExperiment"))
    stop("Input must be a 'TimeSeriesExperiment' object.")
  if (!validObject(object))
    stop("Invalid TimeSeriesExperiment object.")

  clust_params_default <- list(
    dist = "euclidean",
    dynamic = FALSE,
    hclust_params = list(),
    static_cut_params = list(h = 0.5),
    dynamic_cut_params = list()
  )
  clust_params_update  <- clust_params_default
  for (ele in names(clust.params)) {
    if (ele %in% names(clust_params_update)){
      clust_params_update[[ele]] <- clust.params[[ele]]
    }
  }

  group <- cluster <- feature <- NULL
  n.top.feat <- min(n.top.feat, length(rownames(object)))
  if (any(groups.selected == "all")){
    groups.selected <- unique(groups(object))
  }

  if (is.null(timeSeries(object, name = "ts_collapsed"))) {
    object <- collapseReplicates(object)
    object <- makeTimeSeries(object)
  }

  ts_collapsed <- timeSeries(object, "ts_collapsed") %>%
    dplyr::select(-replicate) %>%
    filter(group %in% groups.selected)  # filter to only the chosen groups
  # Find top "n.top.feat" most variable features
  feat_tmps <- ts_collapsed %>%
    dplyr::select(-starts_with("Lag"))
  feat_tmps$sd <- apply(feat_tmps %>% dplyr::select(-feature, -group),
                        1, sd, na.rm = TRUE)

  feat_sd <- feat_tmps %>%
    group_by(group) %>%
    top_n(n = n.top.feat, wt = sd)
  top_features <- unique(feat_sd[["feature"]])

  # Aggregate timecourses across all "groups" and recompute lags
  if(length(groups.selected) > 1) {
    message("Averaging timecourses over all 'groups' selected ",
            "and recomputing lags with coefficients: ",
            paste0(lambda, collapse = " "))
    ts_collapsed <- feat_tmps %>%
      dplyr::select(-group, -sd) %>%
      group_by(feature) %>%
      summarise_all(mean)
    ts_with_lags <- .addLagsToTimeSeries(
      ts_collapsed %>% dplyr::select(-feature), lambda = lambda)
    ts_collapsed <- cbind(feature = ts_collapsed$feature, ts_with_lags)
  } else {
    if(!any(grepl("Lag_", colnames(timeSeries(object, "ts_collapsed"))))){
      object <- addLags(object, lambda = lambda)
    }
    ts_collapsed <- dplyr::select(feat_tmps, -group)
  }
  # Cluster a subset of features
  ts_subset <- ts_collapsed %>%
    filter(feature %in% top_features) %>%
    remove_rownames() %>%
    column_to_rownames("feature")
  clust_params_update[["X"]] <- ts_subset
  res_cluster_subset <- do.call(clusterData, clust_params_update)
  cluster_hclust <- res_cluster_subset$hclust
  clust_map <- res_cluster_subset$clust_map %>%
    dplyr::select(feature, cluster)
  clust_centroids <-  res_cluster_subset$clust_centroids %>%
    remove_rownames() %>%
    column_to_rownames("cluster")

  # Assign the rest of the genes to the closest cluster centroid
  ts_remain <- ts_collapsed %>%
    filter(!feature %in% clust_map$feature) %>%
    remove_rownames() %>%
    column_to_rownames("feature")
  dist_to_nearest_clust <- proxy::dist(ts_remain, clust_centroids)
  clst.remain <- apply(dist_to_nearest_clust, 1, function(x) {
    colnames(dist_to_nearest_clust)[which.min(x)]
  })
  clust_map_remain <- data.frame(
    feature = names(clst.remain),
    cluster = clst.remain)
  clust_map$used_for_hclust <- rep(TRUE, nrow(clust_map))
  clust_map_remain$used_for_hclust <- rep(FALSE, nrow(clust_map_remain))
  final_cluster_map = rbind(clust_map, clust_map_remain)

  freq_df <- final_cluster_map %>%
    group_by(cluster) %>%
    summarise(freq = n()) %>%
    arrange(desc(freq)) %>%
    mutate(cluster_name = paste0("C", seq_len(nrow((.)))))

  final_cluster_map <- final_cluster_map %>%
    mutate(cluster = factor(
      cluster, levels = freq_df$cluster, labels = freq_df$cluster_name))

  # rowData(object) <- suppressMessages(
  #   DataFrame(as.data.frame(rowData(object)) %>%
  #               left_join(final_cluster_map))
  # )

  res_cluster_subset$clust_map <- res_cluster_subset$clust_map %>%
    mutate(cluster = factor(
      cluster, levels = freq_df$cluster, labels = freq_df$cluster_name))

  res_cluster_subset$final_cluster_map <- final_cluster_map
  slot(object, name = "clusterAssignment", check = TRUE) <-
    res_cluster_subset
  return(object)
}


create_c_map_standard_cluster <- function(object,top_feats=500){
  #Here we do the clustering like in the tutorial

  params_for_clustering <- list(
    dynamic = TRUE,
    dynamic_cut_params = list(deepSplit = TRUE, minModuleSize = 20))

  object <- clusterTimeSeries_mod(
    object, n.top.feat = top_feats, groups = "all",
    clust.params = params_for_clustering)
  #Create standard cluster map
  cluster_map <- clusterMap(object)

  hclust_obj <- clusterAssignment(object, "hclust")
  svg("TS_results/TS_clusters_standard_clustering_dendrogramme.svg",width=8,height=5)
  plot(x = hclust_obj, labels = FALSE, xlab = "genes", sub = "")
  dev.off()

  cluster_map<-cluster_map[,1:2]

  plotTimeSeriesClusters(object, transparency = 0.2,scales = 'free_x')+
    scale_color_manual(values = c( 'Control' = "#1f78b4", 'Treated' = "#e31a1c")) +
    theme(strip.text = element_text(size = 10))

  return(cluster_map)
}



conditional_DESeq2_DE <- function(time_object, groups=c('control','experiment'),genes_of_interest=c(),species='hg38'){
  for (item in groups){
    file_sample <- colData(time_object)[colData(time_object)$group==item,]
    file_sample <- file_sample$sample
    file <- assays(time_object)$raw[,colnames(assays(time_object)$raw) %in% file_sample]
    file <-as.data.frame(file)
    file<-add_column(file, Genes = rownames(assays(time_object)$raw), .before = 1)

    if (exists('core_files_1')==T){
      core_files_2 <-file
    }else{
      core_files_1<-file
    }
  }

  gene_num_list <- list()
  for (tp in unique(timepoints(time_object))){
    message("testing timepoint: ", tp)

    #Create name of experiment for summaries
    exp_name<-paste0(groups[1],'_vs_',groups[2],'_TP_',tp)

    target_files<-names(which(timepoints(time_object)==tp))

    t_files_1<-target_files[target_files %in% colnames(core_files_1)]
    t_files_2<-target_files[target_files %in% colnames(core_files_2)]

    tp_file_1<-core_files_1[,c('Genes',t_files_1)]
    tp_file_2<-core_files_2[,c('Genes',t_files_2)]

    gtf_p <- 'data_files/gencode.v35.annotation.gtf'
    #Groups are 'reversed' since control is group 1, and we want the results to be in regards to the experiment
    DE_info<-DESeq2_gtf_wrapper(tp_file_2,tp_file_1,groups[2],groups[1],gtf_path=gtf_p,
                                genes_of_interest=genes_of_interest,
                                result_location = 'TS_results/DE_results_conditional',time_point = tp,species=species)
    gene_num_list[[exp_name]] <- DE_info
  }
  df <- data.frame(matrix(unlist(gene_num_list), nrow=length(gene_num_list), byrow=T))
  rownames(df)=names(gene_num_list)
  colnames(df)=c('Number of Genes', "Number of Significant Genes")
  write.csv(df,"TS_results/DE_results_conditional/DE_summary.csv", row.names = TRUE)
  create_tables_genes_interest(genes_of_interest,names(genes_of_interest),pval_or_padj=p_val_filter_type,
                               path_to_DE = 'TS_results/DE_results_conditional/',results_path = 'TS_results/DE_results_conditional/tables_genes_interest/')
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A wrapper function to plot differential gene expression results as a heatmap
#' The function simply processes the Differential expression results to only send 
#' the correct time points and genes to the plotting function
#' 
#' @param DE_object
#' @param genes_to_plot
#' @param save_name_plot
#' @param save_name_csv
#' @param time_to_plot
#' @param cluster_columns
#' @param save_path
#' @param show_genes_on_plot
#' @param custom_height
#' @param custom_width
#' 
#' 
#' @examples
#' 
#'
#' @export
DE_plot_wrapper_for_heatmap <- function(DE_object,genes_to_plot,save_name_plot,
                                        save_name_csv, time_to_plot =NULL, 
                                        cluster_columns=F, save_path='TS_results',
                                        show_genes_on_plot=F, custom_height=30,
                                        custom_width=25){
  
  #Set up the object for plotting
  sub_tmpDE <- list()
  for (ts in names(DE_object)){
    sub_tmpDE[[ts]] <- DE_object[[ts]][DE_object[[ts]]$feature %in% genes_to_plot,]
  }
  
  #Set up the times that will be plotted
  if (is.null(time_to_plot)==T){
    time_to_plot <- names(sub_tmpDE)
  }
  
  
  svg(paste0(save_path,save_name_plot), height = custom_height, width = custom_width)
  print(plotHeatmap_using_DE_data(main_data=sub_tmpDE,experiment_selection = time_to_plot,
                                  save_name = paste0(save_path,save_name_csv),gene_dend_width=2,
                                  show_gene_symbols = show_genes_on_plot,cluster_the_colums = F))
  dev.off()
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A wrapper function to create enrichment plots for every cluster of the enrich_res
#' results
#' 
#' The functions creates the enrichment dot plot and a heatmap for that cluster,
#' the results are sotred in a 'individual_cluster' folder, then an individual folder
#' per cluster within the 'individual_clusters' folder.
#' 
#' 
#' @param enrich_res The enrichment pathway results from the pathwayEnrichment_custom_cluster_map ()
#' @param cluster_map The cluster map
#'
#' @examples
#' 
#'
#' @export
heat_enrich_wrap_for_clusters <- function(object,enrich_res,c_map,do_heatmap=T,file_name='individual_clusters'){
  
  main_path<-paste0('TS_results/',file_name)
  dir.create(main_path)
  for (clust in names(enrich_res)){
    
    if(nrow(enrich_res[[clust]])==0){
      next
    }
    
    current_path=paste0(main_path,'/',clust)
    dir.create(current_path)
    
    ggsave(paste0(current_path,'/',clust,'.png'),width=13,height=9,
           plotEnrichment(
             enrich = enrich_res[[clust]], n_max = 15) +
             ggtitle(paste0("Cluster , ", clust, " enrichment"))
    )
    if (do_heatmap==T){
      target_exp_plot <- colnames(object)
      genes_from_cluster <- c_map$feature[c_map$cluster==clust]
      plotHeatmap_part_mod(object, num.feat = 0, use_samples= target_exp_plot,
                           specific_genes = genes_from_cluster,
                           heat_name=paste0(current_path,'/variable_heat_',clust),
                           save_name = paste0(current_path,'/variable_genes_',clust,'.csv'))
      
      # temp_var <- TS_object
      # temp_var_2 <- enrich_res
      # load('Result_area/covid/TimeSeries_results/un-unlogged_data_trials/TS_results_all_patients_clin_filter_genes/TS_object.RData')
      # target_exp_plot <- colnames(object)
      # plotHeatmap_part_mod(object, num.feat = 0, use_samples= target_exp_plot,
      #                      specific_genes = genes_from_cluster,
      #                      heat_name=paste0(current_path,'/all_patients_variable_heat_',clust),
      #                      save_name = paste0(current_path,'/all_patients_variable_genes_',clust,'.csv'))
      # 
      # TS_object<-temp_var
      # enrich_res<-temp_var_2
    }
    
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Creates a clustered heat map using PART algorithm from clusterGenomics
#' 
#' clusterGenomics developed a 'deeper' clustering using PART 
#' (Partitioning Algorithm based on Recursive Thresholding). This function illustrates
#' that clustering on a heat map of the top n variable genes in the time series dataset.
#' The heatmap also helps illustrate the clusters through colored bars.
#' 
#' The original data is scaled for a better visualization, then saved into a csv
#' The csv contains the gene symbols, the clusters to which they have been paired
#' as well as the value associated to that gene per sample.
#' 
#' 
#' @param object A time series object
#' @param num.feat The number of features/genes that will be plotted
#' @param scale Weather or not to scale the values
#' @param feat_desc What 'gene_names' to use (feature or symbol)
#' @param sample_desc The name given to the group of samples (by default, sample)
#' @param save_name The save_name to be given to the csv file, put null to not save
#' @param dist_param Parameter for the dist function of clustering, should not be changed
#' if the goal is to preserve part clustering
#' @param hclust_param Parameter for the hclust function of clustering, should not be changed
#' if the goal is to preserve part clustering
#'
#' @examples
#' plotHeatmap_part_mod(TS_object, num.feat = 300, save_name = 'TS_results/top_300_variable_genes.csv')
#'
#' @export
plotHeatmap_part_mod <- function(object, num.feat = 200, scale = T, 
                                 feat_desc = "symbol", sample_desc = "sample",
                                 save_name="plotHeat_part.csv", use_samples=NULL,
                                 specific_genes=NULL, heat_name='my_heat',
                                 dist_param="euclidean", hclust_param="average",
                                 use_seed='123456')
{   
  group <- timepoint <- NULL
  cnts <- assays(object)$norm
  rownames(cnts) <- rowData(object)[, feat_desc]
  colnames(cnts) <- colData(object)[, sample_desc]
  
  found_groups<-setNames(object@group,colnames(object))
  found_replicates<-setNames(object@replicate,colnames(object))
  found_timepoints<-setNames(object@timepoint,colnames(object))
  
  if (is.null(use_samples)==F){
    cnts <- cnts[, use_samples]
    found_groups =found_groups[use_samples]
    found_replicates =found_replicates[use_samples]
    found_timepoints =found_timepoints[use_samples]
  }
  
  if (is.null(specific_genes)==T){
    top_feat <- apply(cnts, 1, sd)
    top_feat <- names(top_feat)[order(-top_feat)[seq_len(num.feat)]]
  }else{
    top_feat <- specific_genes
    gene_not_found <- top_feat[!top_feat %in% rownames(cnts)]
  }
  
  
  top_feat <- top_feat[top_feat %in% rownames(cnts)]
  
  cols_ordered <- order(found_groups, found_replicates, 
                        found_timepoints)
  
  #creates, orders, and scales the normalized count matrix for the top features selected
  Y <- cnts[top_feat, cols_ordered]
  # Y<-Y+1
  # Y<-log10(Y)
  if (scale) {
    Y <- t(scale(t(Y), center = TRUE, scale = TRUE))
  }
  
  #Creates a dataframe to be used as a means to label and illustrate time points,
  #replicates, and groups
  smpdf <- colData(object) %>%
    as.data.frame() %>%
    dplyr::select(group, replicate, timepoint) %>%
    dplyr::arrange(group, replicate, timepoint) %>%
    dplyr::mutate(timepoint = as.numeric(timepoint))
  
  if (is.null(use_samples)==F){
    smpdf <- smpdf[rownames(smpdf) %in% use_samples,]
  }
  
  #sets up the format for the group labeling
  n_group <- length(unique(found_groups))
  cols <- RColorBrewer::brewer.pal(9, name = "Set1")[
    seq_len(min(9, n_group))]
  group_cols <- grDevices::colorRampPalette(colors = cols)(n_group)
  names(group_cols) <- unique(found_groups)
  
  #sets up the format for the replicate labeling
  n_replicates <- length(unique(found_replicates))
  cols <- RColorBrewer::brewer.pal(8, name = "Set3")[
    seq_len(min(8, n_replicates))]
  rep_cols <- grDevices::colorRampPalette(colors = cols)(n_replicates)
  names(rep_cols) <- unique(found_replicates)
  
  #sets up the format for the time point labeling
  time_cols <- circlize::colorRamp2(
    breaks = seq(min(found_timepoints), max(found_timepoints),
                 length.out = 10),
    colors = viridis(10))
  
  
  #Merges the labeling formats made above into a heatmap annotation format
  ha1 <- ComplexHeatmap::HeatmapAnnotation(
    df = smpdf, 
    col = list(group = group_cols, replicate = rep_cols, 
               timepoint = time_cols))
  
  #Calculates the clustering using the 'part' algorithm
  set.seed(use_seed)
  calculated_clusters = part(Y,B=200,minSize=10)
  rowclust = hclust(dist(Y,method=dist_param),method=hclust_param)
  
  
  #Sets up the format for the color bar which will illustrate the different clusters found
  num_clusts <- length(unique(calculated_clusters$lab.hatK))
  cols <- RColorBrewer::brewer.pal(8, name = "Set3")[seq_len(min(8, num_clusts))]
  clust_cols <- grDevices::colorRampPalette(colors = cols)(num_clusts)
  names(clust_cols) <- unique(as.character(calculated_clusters$lab.hatK))
  
  
  # clust_ordered <- as.character(sort.default(unique(calculated_clusters$lab.hatK)))
  clust_ordered <- unique(as.character(calculated_clusters$lab.hatK[rowclust$order]))
  #Cluster illustration stored as rowAnnotation
  row_annot <- ComplexHeatmap::rowAnnotation(gene_cluster = as.character(calculated_clusters$lab.hatK), 
                                             col = list(gene_cluster=clust_cols),
                                             show_annotation_name=F,
                                             annotation_legend_param = list(title = "clusters", at = clust_ordered,
                                                                            labels = clust_ordered))
  
  #Saves the file as a csv if save_name is not null
  if (is.null(save_name)==F){
    FD_2 <- Y[rowclust$order,]
    FD_2 <- as.data.frame(FD_2)
    
    #Need to obtain the features in the same order for the symbols which are?
    # symbols <- as.data.frame(rownames(FD_2))
    # colnames(symbols)<-c('symbol')
    # TS_obj_feat_symbol_data<- (as.data.frame(rowData(TS_object)))
    # features <- TS_obj_feat_symbol_data[TS_obj_feat_symbol_data$symbol %in% symbols$symbol,]
    # features <- features[match(symbols$symbol,features$symbol),]
    FD_2 <- add_column(FD_2, gene_cluster = calculated_clusters$lab.hatK[rowclust$order], .before = 1)
    # FD_2 <- add_column(FD_2, features = features$feature, .before = 1)
    write.csv(FD_2,save_name, row.names = T)
  }
  
  
  svg(paste0(heat_name,"_with_names.svg"),height=30,width=20)
  print(ComplexHeatmap::Heatmap(
    Y, name = "Z-score", cluster_columns = FALSE,cluster_rows=rowclust,
    top_annotation = ha1, row_names_gp = grid::gpar(fontsize = 8),left_annotation = row_annot,
    show_row_names = T))
  dev.off()
  
  
  
  svg(paste0(heat_name,".svg"))
  print(ComplexHeatmap::Heatmap(
    Y, name = "Z-score", cluster_columns = FALSE,cluster_rows=rowclust,
    top_annotation = ha1, row_names_gp = grid::gpar(fontsize = 8),left_annotation = row_annot,
    show_row_names = F))
  dev.off()
  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Creates a clustered heat map using PART algorithm from clusterGenomics for 
#' differential expression data from the TimeSeriesExperiment package
#' 
#' clusterGenomics developed a 'deeper' clustering using PART 
#' (Partitioning Algorithm based on Recursive Thresholding). This function illustrates
#' that clustering on a heat map of the top n variable genes in the time series dataset.
#' The heatmap also helps illustrate the clusters through colored bars.
#' 
#' The heatmap, by default, plots the log fold change between the time series.
#' The function takes in the necessary time series (as inputted by the user) and 
#' subsets the data to only include those time series, the genes are then balanced so
#' that only genes appearing in the selected time points are plotted.
#' 
#' The function can also save the clustering results along with the values
#' in the format of a csv file
#' 
#' 
#' @param main_data The differential expression matrix to be plotted
#' @param experiment_selection A vector of names representing the time points to be plotted
#' @param data_to_use The data type that will be plotted
#' @param dist_param Parameter for the dist function of clustering, should not be changed
#' if the goal is to preserve part clustering
#' @param hclust_param Parameter for the hclust function of clustering, should not be changed
#' if the goal is to preserve part clustering
#' @param save_name The name of the csv file containing the clustering and heat map date
#' (set to NULL to not save the csv)
#' @param gene_dend_width The width of the row side dendrogram
#' @param show_gene_symbols A boolean to check if gene symbols should be illustrated on the right
#' hand side of the heat map
#'
#' @examples
#' exp_select <- c('2.5','4','6','9','13')
#' png("TS_results/DE_heat_map.png", 750, 500)
#' plotHeatmap_using_DE_data(main_data=tmp_de,experiment_selection = exp_select,
#'                           save_name = 'TS_results/DE_clust_heat.csv',gene_dend_width=2,
#'                           show_gene_symbols = F)
#' dev.off()
#'
#' @export
plotHeatmap_using_DE_data <- function(main_data,experiment_selection,data_to_use='logFC',
                                      dist_param="euclidean",hclust_param="average",
                                      save_name="DE_heat.csv",gene_dend_width=2,
                                      show_gene_symbols=T,cluster_the_colums=F, ...)
{   
  #Sets the columns to keep
  cols_needed=c('symbol','logFC')
  
  #subsets the main data to only contain the necessary time points and genes
  for (item in names(main_data)){
    if (item %in% experiment_selection){
      sub_data <- main_data[[item]]
      sub_data <- sub_data[cols_needed]
      names(sub_data)[names(sub_data) == 'logFC'] <- item
      if (exists("final_df")==F){
        final_df <- sub_data
      }else{
        final_df <- merge(final_df,sub_data,by='symbol')
      }
    }
    
  }
  rownames(final_df)<-final_df$symbol
  final_df <- final_df[,-1]
  
  #Convert to numerical
  final_df=as.matrix(final_df)
  calculated_clusters = part(final_df,B=100)
  
  #cluster rows and columns using 'part'
  rowclust = hclust(dist(final_df,method=dist_param),method=hclust_param)
  colclust = hclust(dist(t(final_df), method=dist_param),method=hclust_param)
  
  #Creates an annotation bar to help visualize the clusters
  num_clusts <- length(unique(calculated_clusters$lab.hatK))
  cols <- RColorBrewer::brewer.pal(8, name = "Set3")[seq_len(min(8, num_clusts))]
  clust_cols <- grDevices::colorRampPalette(colors = cols)(num_clusts)
  names(clust_cols) <- unique(as.character(calculated_clusters$lab.hatK))
  
  #Sets the cluster annotation bar as a rowAnnotation object
  row_annot <- ComplexHeatmap::rowAnnotation(clusters = as.character(calculated_clusters$lab.hatK), 
                                             col = list(gene_cluster=clust_cols),
                                             show_annotation_name=F)
  
  #Creates the annotation bar for the time points
  num_series <- length(unique(experiment_selection))
  cols <- RColorBrewer::brewer.pal(8, name = "Set1")[seq_len(min(8, num_series))]
  series_col <- grDevices::colorRampPalette(colors = cols)(num_series)
  names(series_col) <- unique(as.character(experiment_selection))
  #Stores the time point annotation as HeatmapAnnotation
  
  col_annot <- ComplexHeatmap::HeatmapAnnotation(exp_series = as.character(experiment_selection), 
                                                 col = list(exp_series=series_col),
                                                 show_annotation_name=FALSE,show_legend = T,
                                                 annotation_legend_param = list(title = "timepoints", at = names(series_col), 
                                                                                labels = names(series_col)))
  
  #If asked, saves the information in a csv file
  if (is.null(save_name)==F){
    # Order data matrix according to order in clustering og plot heatmap:
    if (cluster_the_colums==T){
      FD_2 <- final_df[rowclust$order, colclust$order]
    }else{
      FD_2 <- final_df[rowclust$order, ]
    }
    
    FD_2 <- final_df[rowclust$order,]
    FD_2 <- as.data.frame(FD_2)
    FD_2 <- add_column(FD_2, gene_cluster = calculated_clusters$lab.hatK[rowclust$order], .before = 1)
    write.csv(FD_2,save_name, row.names = T)
  }
  
  if (cluster_the_colums ==T){
    ComplexHeatmap::Heatmap(
      final_df,name=data_to_use, cluster_columns = colclust,
      row_names_gp = grid::gpar(fontsize = 8),cluster_rows=rowclust,
      row_dend_width = unit(gene_dend_width, "cm"),left_annotation = row_annot,
      top_annotation = col_annot,show_column_names = F,show_row_names = show_gene_symbols)
  }else{
    ComplexHeatmap::Heatmap(
      final_df,name=data_to_use,
      row_names_gp = grid::gpar(fontsize = 8),cluster_rows=rowclust,cluster_columns = F,
      row_dend_width = unit(gene_dend_width, "cm"),left_annotation = row_annot,
      top_annotation = col_annot,show_column_names = F,show_row_names = show_gene_symbols)
  }
  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' 
#' 
#' 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#'
#' @examples
#' 
#'
#' @export
.addLagsToTimeSeries <- function(timeseries, lambda) {
  if(is.null(lambda)) stop("Need to specify weights for lags.")
  nT <- ncol(timeseries)
  time_names <- colnames(timeseries)
  timeseries <- as.matrix(timeseries)
  lags <- lapply(seq_along(lambda), function(i) {
    ilag <- lambda[i] * t(diff(t(timeseries), lag = i))
    colnames(ilag) <- paste0(
      "Lag_", time_names[seq((i+1), nT)], "_", time_names[seq(1,(nT-i))])
    return(ilag)
  })
  lags <- do.call("cbind", lags)
  res <- as.data.frame(cbind(timeseries, lags))
  return(res)
}




library(edgeR)
library(limma)
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' 
#' 
#' 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#'
#' @examples
#' 
#'
#' @export
timepointDE_mod <- function(object, timepoints = "all", 
                            min_gene_sum = 1, alpha = 0.05) 
{
  if (!is(object, "TimeSeriesExperiment")) 
    stop("Input must be a 'TimeSeriesExperiment' object.")
  if (!validObject(object)) 
    stop("Invalid TimeSeriesExperiment object.")
  feature <- group <- time <- NULL
  if (any(timepoints == "all")) {
    timepoints <- sort(unique(timepoints(object)))
  }
  if(!all(timepoints %in% unique(timepoints(object)))) {
    stop("One or more entries of 'timepoints' not found in ",
         "timepoints(object).")
  }
  # we rename "group" to "condition" beacuse in 'model.matrix()' misuses
  # 'group' in design argument
  sample_data <- colData(object) %>%
    as.data.frame() %>%
    rename(condition = group)
  feature_data <- rowData(object) %>%
    as.data.frame()
  
  limma.res <- lapply(timepoints, function(t) {
    message("testing timepoint: ", t)
    smp.t <- sample_data %>% filter(time == t)
    x.t <- assays(object)$raw[, smp.t$sample]
    # filter out very sparse genes:
    x.t <- x.t[rowSums(x.t) >= min_gene_sum, ]  
    gdata.t <- feature_data %>% filter(feature %in% rownames(x.t))
    
    # Construct a DGEList object
    dge <- DGEList(counts = x.t, genes = gdata.t, samples = smp.t)
    
    # Compute size factor (lib sizes)
    dge <- calcNormFactors(dge)
    
    # Construct a design model matrix
    dsgn <- model.matrix(~ condition, dge$samples)
    
    # Voom normalization
    v <- voom(counts = dge, design = dsgn, plot = FALSE)
    
    # Fit the model
    fit <- lmFit(v, dsgn)
    
    # Compute test statistics
    fit <- eBayes(fit)
    
    # Statistics table
    t.res.limma <- suppressMessages(
      topTable(fit, adjust.method="BH", number = Inf, p.value = alpha)
    )
    return(t.res.limma)
  })
  names(limma.res) <- timepoints
  
  diff_exp <- differentialExpression(object)
  diff_exp$timepoint_de <- limma.res
  slot(object, name = "differentialExpression", check = TRUE) <- diff_exp
  return(object)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' 
#' 
#' 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#'
#' @examples
#' 
#'
#' @export
my_diff_exp_TS_func <- function(object,conditions){
  main_dir<-getwd()
  setwd('TS_results')
  
  TP_vect<-object@timepoint
  names(TP_vect)=colnames(object)
  
  GP_vect<-object@group
  names(GP_vect)=colnames(object)
  
  
  
  for (tp in unique(TP_vect)){
    TP_subset<-TP_vect[TP_vect==tp]
    GP_subset<-GP_vect[names(GP_vect) %in% names(TP_subset)]
    
    file_control<-names(GP_subset)[GP_subset==group_names[1]]
    file_control <- assays(time_object)$raw[,colnames(assays(object)$raw) %in% file_control]
    file_control <-as.data.frame(file_control)
    file_control<-add_column(file_control, Genes = rownames(assays(object)$raw), .before = 1)
    
    file_experiment<-names(GP_subset)[GP_subset==group_names[2]]
    file_experiment <- assays(time_object)$raw[,colnames(assays(object)$raw) %in% file_experiment]
    file_experiment <-as.data.frame(file_experiment)
    file_experiment<-add_column(file_experiment, Genes = rownames(assays(object)$raw), .before = 1)
    
    
    View(file_control)
    storage<-DESeq2_gtf_wrapper(File_1=file_control, File_2=file_experiment,
                                cond_1=group_names[1], cond_2=group_names[2], gtf_path=NULL,
                                genes_of_interest=c('AICDA'),species='Mm')
    print(storage)
    break
  }
  
  
  
  setwd(main_dir)
  
  
}
