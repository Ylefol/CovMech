library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)

# Process initial data -------------------


#' This function merges all the histological, sample, and patient data into a single file
#' 
#' The function is tailored specifically to the files given for the CRC dataset
#' and should not be used for other sets of clinical or histological files.
#' 
#' @param histo_path file path to the histological data
#' @param barcode_path file path to the barcode information
#' @param patient_path file path to the patient data
#' 
#' @return master_file The merged information from the three files provided
#'
#' @examples
#'
#' @export
create_master_patient_sample_data <-function(histo_path,barcode_path,patient_path){
  hist_data <- read.csv(histo_path)
  hist_data <- hist_data[hist_data$inflammation_status != 9,]#Remove unspecified inflammation status
  
  barcode_info <-read.csv(barcode_path)
  #Match the formatting between the two files
  hist_data$study_ID<-gsub(pattern = ' ',replacement = '',x = hist_data$study_ID)
  hist_data$study_ID<-gsub(pattern = '-',replacement = '',x = hist_data$study_ID)
  barcode_info$Biopsi<-gsub(pattern = 'D ',replacement = 'D',x = barcode_info$Biopsi)
  barcode_info$Biopsi<-gsub(pattern = ' ',replacement = '-',x = barcode_info$Biopsi)
  
  
  #Create the barcode+position (SampleID) and the custom ID which will serve as a means of merging the two files
  barcode_info$SampleID <- paste(barcode_info$Barcode, barcode_info$Position, sep='_')
  barcode_info$customID <- paste(barcode_info$Patient, barcode_info$Biopsi, sep='_')
  hist_data$customID <- paste(hist_data$study_ID,hist_data$nr_FreshFrozen,sep='_')
  
  #Remove any spaces in customIDs
  barcode_info$customID<-gsub(pattern = ' ',replacement = '',x = barcode_info$customID)
  hist_data$customID<-gsub(pattern = ' ',replacement = '',x = hist_data$customID)
  
  
  #Keep only columns of interest
  barcode_info_subset <- barcode_info[,c('SampleID','customID','Location')]
  temp<-barcode_info_subset
  
  hist_data_subset <- hist_data[,c('study_ID','progressor_status','segment','inflammation_status',
                                   'chronic_grade','active_grade','neoplasia_status',
                                   'atrofi_status','kryptDist_status','customID')]
  #Perform the merge, include files which do not have expression data (SampleID)
  hist_and_barcodes <- merge(barcode_info_subset,hist_data_subset,by='customID',all=TRUE)
  
  
  #Prep patient_info for merger
  patient_info <-read.csv(patient_path)
  patient_info$study_ID <- gsub(pattern = ' ',replacement = '',x = patient_info$study_ID)
  cols_to_keep_patient <- c('study_ID','neoplasia_peak','sex','age_at_diagnosis','age_at_inclusion',
                            'UC_CD_IBDU','UC_type','duration_yr','CRC')
  patient_info <- patient_info[,cols_to_keep_patient]

  master_file <- merge(hist_and_barcodes,patient_info,by='study_ID')

  #Use the Status column as a numeric for progressor status
  master_file$Status <- gsub(pattern = 'NP',replacement = 0,master_file$progressor_status)
  master_file$Status <- gsub(pattern = 'P',replacement = 1,master_file$progressor_status)
  
  ###REVIEW THE ORDERING
  new_order <- c(1,2,3,6,5,4,7,15,16,17,20,21,14,18,19,8,9,10,11,12,13)
  master_file <- master_file[,new_order]
  
  
  #Some specific tweaks. There are two cells in the Location column with non-standard
  #labeling. The below code standardizes them
  master_file[master_file=="Sigmoid (Rectum)"]<-'Rectum'
  master_file[master_file=="Illium/Ascending"]<-'Ascending'
  
  return(master_file)
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' A function which reads and pre-processes microarray data
#' 
#' The function will read all microarray files in a given folder, the control probes
#' will be removed. A background correct is performed followed by normalization.
#' 
#' @param path_microarr_data The path to the folder containing the microarray files
#' 
#' @return rna_biop_data The Elist object with the normalized microarray data
#'
#' @examples
#'
#' @export
read_raw_data_CRC <- function(path_microarr_data){
  raw_files0 <- list.files(path_microarr_data, pattern='^AHUS') # list of file names
  
  # Don't include this one
  raw_files <- raw_files0[!raw_files0 == 'AHUS_257236327373_S01_GE1_107_Sep09_2_2.txt']
  
  raw_files_path <- paste(path_microarr_data, raw_files, sep='/') # add path to filenames
  length(raw_files_path)
  
  # load data (see section 4.5 in the limma user guide).
  rna_rawdata <- read.maimages(file=raw_files_path, source='agilent', green.only=TRUE)
  
  #The above creates and ElistRaw object, but what is in here?
  # So the main component is the numeric matrix, called E. It will contain the expression
  # values. In and EListRaw object these expression are unlogged while in a EList object they are log2 values
  # rows correspond to probes and columns to samples
  
  #there is also an Eb object, which contains the unlogged background expression values of the same dimensions as E.
  
  
  # This is what is done in the babelomics suite. 
  #This finds all the control probes, i-e probes labelled as 1 in the controlType columnb
  ctrl_idx <- rna_rawdata@.Data[[4]]$ControlType != 0
  
  #We then remove the control rows from the EListRaw object
  rna_rawdata_noctrl <- rna_rawdata[!ctrl_idx,]
  
  # correct and normalize
  #Concretely, what does this do?
  # This will use the Eb component of the Elist object to correct the E component
  # Recall that the Eb is the unlogged background expression while the E is the unlogged
  # gene expression
  
  # We then normalize between the arrays so that intensities or log ratios have
  # similar distribution across a set of arrays.
  
  data_bckcor <- backgroundCorrect(rna_rawdata_noctrl, method="normexp", offset=16)
  rna_biop_data <- normalizeBetweenArrays(data_bckcor, method="quantile") 
  return(rna_biop_data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A function which renames each microarray column in the Elist with a SampleID
#' 
#' The function processes the file name to extract the SampleID e.g. 257236323159_1_3
#' where the first set of digits is the barcode and the two last number preceded by
#' underscores are the location within the array. Giving this type of ID creates
#' a unique ID per experiment.
#' 
#' This ID is also present in the master file, the function ensures that the names
#' given in the Elist match the names from the master file
#' 
#' @param data The Elist object
#' @param sampledta The dataframe containing the SampleIDs to check with
#' 
#' @return data the updated Elist object
#'
#' @examples
#'
#' @export
give_sample_IDs <- function(data,sampledta){
  sampID <-c()
  for (file_name in colnames(data)){
    split_var <- strsplit(file_name,'_')[[1]]
    new_sampID <- paste0(split_var[6],'_',split_var[11],'_',split_var[12])
    sampID <- c(sampID,new_sampID)
  }
  sampledta <- na.omit(sampledta)

  if (length(sampID[!sampID %in% sampledta$SampleID])>0){
    print("The names do not match with the sample file")
  }else{
    colnames(data$E) <- sampID
    return(data)

  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' barcode removal function
#' 
#' simple function to remove a specific barcode from the Elist and the patient data
#' 
#' @param microarr_data The E list object to be updated
#' @param patient_data Dataframe of the patient data
#' @param barcode_to_remove a string of the barcode to be removed 
#' 
#' @return return_list a list containing the updated Elist and patient_data
#'
#' @examples
#'
#' @export
remove_specific_barcode <- function(microarr_data,patient_data,barcode_to_remove){
  # barcode_to_remove <- '257236327378_1_2'
  patient_data %>% 
    filter(SampleID != barcode_to_remove) -> patient_data
  
  microarr_data <- microarr_data[,!colnames(microarr_data$E) == barcode_to_remove]
  
  return_list <- list(microarr_data,patient_data)
}



prep_matrix <- function(contrast,file_location,sample_dta){
  file_names <- gsub(pattern = '\\.',replacement = '_',x = strsplit(contrast,'_vs_')[[1]])

  file_1<-read.csv(paste0(file_location,file_names[1],'.csv'))
  file_2<-read.csv(paste0(file_location,file_names[2],'.csv'))
  
  all_samples<-c(file_1$SampleID,file_2$SampleID)
  
  g1_vect<-gsub(pattern = F,replacement = 0,x = sample_dta$SampleID %in% file_1$SampleID)
  g1_vect<-gsub(pattern = T,replacement = 1,x = g1_vect)
  

  g2_vect<-gsub(pattern = F,replacement = 0,x = sample_dta$SampleID %in% file_2$SampleID)
  g2_vect<-gsub(pattern = T,replacement = 1,x = g2_vect)
  
  matrix_df<-data.frame(Group1=g1_vect, Group2=g2_vect)
  
  mm <- model.matrix(~Group1+Group2,matrix_df)
  colnames(mm)=c('X.Intercept','G1','G2')
  
  return(list(mm,all_samples))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Wrapper function for a standard differential gene expression run using limma
#' 
#' @param limma_object An Elist object containing the normalized and corrected microarray data
#' @param comparison_list The list of comparisons to run through the limme DE pipeline
#' Comparisons are formatted like so: exp.cov_vs_control.cov
#' Names within the comparison cannot contain underscores
#' @param subset_path The file path to the location of the subsets used in the limma analysis
#' @param sample_dta A file containing SampleIDs which match the order in the limma object
#' @param gene_annot Gene annotation file if gene annotation is needed- can be NULL
#' @param gtf_path path to gtf file used to add chromosomal location information - can be null
#' @param filter_choice either padj or pvalue, serves as the value to filter for significance
#' @param genes_interest a vector contianing the gene symbols of interest (if any)
#' @param species the code for the species analyzed--currently supports hg38,mm10,ce11
#' 
#' @return None
#'
#' @examples
#'
#' @export
limma_DE_wrapper <-function(limma_object,comparison_list,subset_path,sample_dta,gene_annot=NULL,gtf_path=NULL,
                            filter_choice='padj',genes_interest=c(),species='hg38'){
  
  #Create folder for results and checkpoint storage
  dir.create('DE_results')
  dir.create('DE_results/Robjects_checkpoints')
  
  #Create list for save results
  eb_res_list<-list()
  #Create list to store then print the number of significant genes in each comparison
  gene_num_list<-list()
  for (item in comparison_list){
    #Create matrix 
    return_list<-prep_matrix(item,subset_path,sample_dta)
    mm<-return_list[[1]]
    all_samp_IDs<-return_list[[2]]
    
    #Calculate Empirical Bayesian statistic
    eb_res<-calculate_EB(limma_object,mm,item)

    #Create DESeq2 style results, store total number of genes and significant genes
    temp_gene_num_lst <- produce_DE_results_for_limma(eb_res=eb_res,limma_obj = limma_object,used_samples = all_samp_IDs,gene_annot_file=gene_annot,gtf_path=gtf_path,
                                                      filter_choice=filter_choice,gene_int=genes_interest,
                                                      species=species)
    
    #Store number of total and significant genes in list
    gene_num_list[[names(temp_gene_num_lst)]]<-temp_gene_num_lst[[names(temp_gene_num_lst)]]
    
    #annotate/get gene symbols
    if (is.null(gene_annotation)==F){
      eb_res$genes<-perform_annotation(eb_res$genes,gene_annotation)
    }else{
      eb_res$genes$gene_id <-eb_res$genes$GeneName
    }
    
    #Supplement missing gene_ids with systemic names
    eb_res$genes$gene_id <- ifelse(eb_res$genes$gene_id=='', eb_res$genes$SystematicName, eb_res$genes$gene_id)
    #Store eb_res results in a list, will be used for fgsea analysis.
    eb_res_list[[item]]<-eb_res
  }
  #Save all empirical bayesian results
  save(eb_res_list, file='DE_results/Robjects_checkpoints/final_results.rdata', version = 2)
  
  #Create a summary table showing total number of DE genes and number of significant DE genes
  create_DE_summary_table(gene_num_list)
  
  #Creates the results for the genes of interest
  if (length(genes_interest)>0){
    if (is.null(gtf_path)==T){
      target_file<-'DE_raw_data.csv'
    }else{
      target_file<-'DE_results_with_gtf.csv'
    }
    
    create_tables_genes_interest(genes_interest,comparison_list,pval_or_padj=filter_choice,DE_name = target_file)
  }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Function which uses a model matrix to calculate the empirical Bayesian results
#' from a provided microarray normalized and corrected Elist.
#' 
#' @param micro_arr_dta An Elist object containing the normalized and corrected microarray data
#' @param matrix_model the matrix model as created by the \code{prep_matrix()}
#' @param comparison The name of the comparison being run
#' 
#' @return eb_res The empirical Bayesian results
#'
#' @examples
#'
#' @export
calculate_EB <- function(micro_arr_dta,matrix_model,comparison){
  if(is.fullrank(matrix_model)==F){
    matrix_model<-matrix_model[,1:2]
    #Perform linear regression and set-up contrast
    limma_res <- lmFit(micro_arr_dta, matrix_model)
    #Results will be given in accordance to the second element, we want G2 to be the 'control/reference'
    ctrst <- makeContrasts(G1, levels = colnames(coef(limma_res)))
  }else{
    #Perform linear regression and set-up contrast
    limma_res <- lmFit(micro_arr_dta, matrix_model)
    #Results will be given in accordance to the second element, we want G2 to be the 'control/reference'
    ctrst <- makeContrasts(G1 - G2, levels = colnames(coef(limma_res)))
  }
  
  colnames(ctrst) <- c(comparison)
  #Compute Contrasts From Linear Model Fit
  ctrst_res <- contrasts.fit(limma_res, ctrst)
  #Empirical Bayes Statistics For Differential Expression
  eb_res <- eBayes(ctrst_res)
  return (eb_res)
}


#' Function which performs the differential gene expression analysis using empirical 
#' bayesian results as obtained from limma. the function then plots volcano and MA 
#' plots. If a gtf path is provided, the function will also filter the results
#' for genes present in the gtf file, it will then plot the karyoplots for those 
#' same genes.
#' 
#' @param eb_res The empirical Bayesian results
#' @param gene_annot_file The gene annotation file - can be NULL
#' @param gtf_path The path to the gtf file - can be NULL
#' @param L2FC_thresh The log2 fold change threshold for significance
#' @param p_thresh The padj or pvalue threshold used
#' @param filter_choice Either pvalue or padj, indicates the p value type to filter on
#' @param save_path The path in which the results will be saved
#' @param gene_int A vector containing the gene symbols for the genes of itnerest (if any)
#' @param species The code for the species being analyzed. This is only relevant if a gtf file is provided
#' 
#' @return eb_res The empirical Bayesian results
#'
#' @examples
#'
#' @export
produce_DE_results_for_limma <- function(eb_res,limma_obj,used_samples,gene_annot_file=NULL,gtf_path=NULL,L2FC_thresh=1,
                                         p_thresh=0.05,filter_choice='padj',save_path='DE_results/',
                                         gene_int=c(),species='hg38'){
  
  sig_gene_summary<-list()
  for (ii in 1:ncol(eb_res)){
    fname <- colnames(eb_res)[(ii)]
    dir.create(paste0(save_path,fname))
    exp_save_path<-paste0(save_path,fname)
    
    
    eb_res %>%
      topTable(coef=ii, number = Inf)  -> DE_raw
    
    #Add a gene ID column using annotation file provided by agilent
    if (is.null(gene_annot_file)==F){
      DE_raw<-perform_annotation(DE_raw,gene_annot_file)
    }else{
      DE_raw %>% 
        dplyr::rename(
          gene_id = GeneName,
        ) -> DE_raw
    }
    
    
    #Supplement missing gene_ids with systemic names
    DE_raw$gene_id <- ifelse(DE_raw$gene_id=='', DE_raw$SystematicName, DE_raw$gene_id)
    #Filter for columns that we want
    DE_raw <- DE_raw[,c('ProbeName','gene_id','logFC','AveExpr','P.Value','adj.P.Val')]
    
    #Change column names to match DE results that have been created before to allow 
    #these results to be sent to the rest of the pipeline
    colnames(DE_raw)=c('ProbeName','gene_id','log2FoldChange','baseMean','pvalue','padj')
    
    DE_raw =  DE_raw[order(DE_raw$gene_id,DE_raw$padj),]
    DE_raw = DE_raw[ !duplicated(DE_raw$gene_id), ]
    
    
    #Add patients used to the file
    sample_counts<-limma_obj$E[,used_samples]
    genes_for_sort<-limma_obj$genes[,c('ProbeName','SystematicName')]
    row.names(genes_for_sort)=genes_for_sort$ProbeName
    genes_for_sort<-genes_for_sort[DE_raw$ProbeName,]
    sample_counts<-sample_counts[genes_for_sort$SystematicName,]
    
    DE_raw<-cbind(DE_raw,sample_counts)
    #Save the DE_raw results
    write.csv(DE_raw, file=paste0(exp_save_path,'/DE_raw_data.csv'), row.names=FALSE,quote = F)
    
    #Create gtf-like results
    pval_for_naming<-paste(unlist(str_split(p_thresh,'\\.')),collapse='')
    if (is.null(gtf_path)==F){
      DE_raw <- add_gtf_to_DE_res(gtf_path, DE_results=DE_raw,folder_n=exp_save_path)
      sig_file_name<-paste0('DE_gtf_significant_',L2FC_thresh,'_',pval_for_naming,'.csv')
    }else{
      sig_file_name<-paste0('DE_significant_',L2FC_thresh,'_',pval_for_naming,'.csv')
    }
    
    total_gene_number <- nrow(DE_raw)
    ggsave(paste0(exp_save_path,"/volcano_plot_with_top_10.png"),dpi=300,width=21, height=19, units='cm',
           volcanoplot_alt(DE_res = DE_raw,filter_choice = filter_choice,
                           plot_title = fname,label_top_n=10,FDR_thresh = p_thresh, l2FC_thresh = L2FC_thresh,genes_of_interest = gene_int)
    )
    ggsave(paste0(exp_save_path,"/volcano_plot.png"),dpi=300,width=18, height=19, units='cm',
           volcanoplot_alt(DE_res = DE_raw,filter_choice = filter_choice,
                           plot_title = fname,label_top_n=0,FDR_thresh = p_thresh, l2FC_thresh = L2FC_thresh,genes_of_interest = gene_int)
    )
    ggsave(paste0(exp_save_path,"/ma_plot.png"),dpi=300,width=21, height=19, units='cm',
           plot=maplot_alt(DE_res = DE_raw,filter_choice = filter_choice,
                           plot_title = fname,FDR_thresh = p_thresh, l2FC_thresh = L2FC_thresh,genes_of_interest = gene_int)
    )

    sig_raw <- filter_for_significance(DE_raw,l2fc_thresh = L2FC_thresh, 
                                       filter_thresh = p_thresh,filter_choice=filter_choice,
                                       name_new_file = paste0(exp_save_path,"/",sig_file_name))
    
    #Records the number of genes in the file, i.e the number of significant genes in our dataset
    sig_gene_number <- nrow(sig_raw)
    
    #Stores the total number of genes and the number of significant genes in a list with the experiments name
    sig_gene_summary[[fname]]<-c(total_gene_number,sig_gene_number)
    
    if (is.null(gtf_path)==F){
      #Set up for karyoplotter
      if (species=='hg38'){
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      }else if (species=='mm10'){
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      }else if(species=='ce11'){
        library(TxDb.Celegans.UCSC.ce11.refGene)
        txdb <- TxDb.Celegans.UCSC.ce11.refGene
      }
      all.genes <- genes(txdb)
      
      # A try catch loop to generate all the karyoplots
      # We use ggsave to save the function, this will always cause a grid.draw error
      # The grid.draw error will not affect the outcome of the plot, but will halt the code
      # Hence the try catch
      karyo_folder_name=paste0(exp_save_path,"/karyoplots")
      dir.create(karyo_folder_name)
      for (chr in unique(DE_raw$seqnames)){
        svg(paste0(karyo_folder_name,"/",chr,".svg"),width=8.5)
        plot_DE_karyotype(DE_raw,genome_to_plot = species, chromo_submitted = c(chr), filter_choice = filter_choice,
                          known_genes=all.genes,genes_of_interest = gene_int)
        dev.off()
      }
    }


    save(DE_raw,file=paste0(save_path,"Robjects_checkpoints/deseq2_",fname,".RData"))#Saves all objects in memory in a binary file
    
  }
  return(sig_gene_summary)
}



#' Function which creates a summary table for the set of comparisons performed.
#' The table contains the total number of genes in each comparison as well as the number
#' of significant genes in each comparison,
#' The function also warns the user if two comparison have the same number of significant
#' DEGs as this may indicated comparison duplication.
#' 
#' @param gene_num_list List containing the total number of genes for each comparison
#' done, along with the number of significant genes in each comparison
#' 
#' @return eb_res The empirical Bayesian results
#'
#' @examples
#'
#' @export
create_DE_summary_table <- function(gene_num_list,save_path='DE_results/'){
  #Save the total number of genes and number of significant genes per analysis to a csv
  df <- data.frame(matrix(unlist(gene_num_list), nrow=length(gene_num_list), byrow=T))
  rownames(df)=names(gene_num_list)
  colnames(df)=c('Number of Genes', "Number of Significant Genes")
  write.csv(df,paste0(save_path,"DE_summary.csv"), row.names = TRUE)
  
  if (length(unique(df$`Number of Significant Genes`)) != length(df$`Number of Significant Genes`)){
    message("Warning: Two or more of the results have the same number of significant DEGs")
  }
}


filter_data_based_on_subsets <-function(microarr_data,target_file,comp_list,subset_path){
  subsets<-c()
  for (i in comp_list){
    for (sub in strsplit(i,'_vs_')[[1]]){
      subsets<-c(subsets,sub)
    }
  }
  subsets<-unique(subsets)
  subsets<-gsub(x = subsets,pattern = '\\.',replacement = '_')
  print(subsets)
  
  IDs_in_subsets<-c()
  for (subset in subsets){
    temp_file<-read_csv(paste0(subset_path,subset,'.csv'))
    IDs_in_subsets<-c(IDs_in_subsets,temp_file$Array)
  }
  IDs_in_subsets<-unique(IDs_in_subsets)
  print(length(IDs_in_subsets))
  
  target_file<-target_file[target_file$Array %in% IDs_in_subsets,]
  IDs_in_target<-target_file$ID
  microarr_data<-microarr_data[,colnames(microarr_data$E) %in% IDs_in_target]
  
  return(list(target_file,microarr_data))
}

#' Retain and save significant DE as per user specifications
#'
#' Filter the given differential expression analysis results using log2FoldChange 
#' and p-value. The absolute value of the log2Foldchange must be greater
#' or equal to the given threhold to pass the filter. The padj (adjusted p-value) 
#' must be lower than the given p-value threshold to pass the filter.
#' 
#' The function then saves the result to a csv file which is named according to user specifications
#'
#' @param DE_file The path to the differential expression file
#' @param l2fc_thresh The threshold for log2FoldChange
#' @param pval_thresh The threshold for padj (adjusted-pvalue)
#' @param name_new_file The name which will be given to the generated results. **Name must contain '.csv' at the end**
#' 
#' @return None
#'
#' @examples
#' filter_for_significance(DE_file, 2,0.05,'DE_gtf_significant.csv')
#'
#' @export
filter_for_significance <- function (DE_file, l2fc_thresh, filter_thresh,filter_choice, name_new_file){
  #Reads the DE file given
  # DE_file <-read.csv(file = DE_file,check.names=FALSE)
  
  #Filters based on the log2FoldChange using absolute value to preserve down- and up-regulated significant genes
  DE_file <- DE_file[abs(DE_file$log2FoldChange) >= l2fc_thresh,]
  #Filters based on the p-value, padj (adjusted p-value) must be lower than the given threshold to be kept
  DE_file <- DE_file[DE_file[filter_choice] < filter_thresh,]
  
  #Saves the file to a csv as per the users specifications
  write.csv(DE_file, file=name_new_file,row.names=FALSE,quote = FALSE)
  
  return(DE_file)
  
}

my_PCA_plotting_function_custom<-function(filtered_array_dat,color_vector,label_vect,folder_name='my_PCAs'){
  dir.create(folder_name)
  filtered_array_dat$E<-log2(filtered_array_dat$E)
  pdf(file=paste0(folder_name,'/OG_health.pdf'))#Prevent plotting of the MDS, just want the values
  arraydata.pca <- plotMDS(filtered_array_dat, pch = 16, gene.selection = "common", col=color_vector)
  dev.off()
  PC1 <- paste("PC1 (", round(arraydata.pca$var.explained[1]*100, 1), "%)", sep = "")
  PC2 <- paste("PC2 (", round(arraydata.pca$var.explained[2]*100, 1), "%)", sep = "")
  percentVar<-c(PC1,PC2)
  
  pca_data<-data.frame(PC1=arraydata.pca$x,PC2=arraydata.pca$y,
                       group=factor(names(label_vect),levels=unique(names(label_vect))),
                       colors=factor(unname(label_vect),levels=unique(unname(label_vect))))
  
  
  # plotMDS(filtered_array_dat, pch = 16, gene.selection = "common", col=color_vector)
  # legend(par("usr")[2], par("usr")[4], unique(names(label_vect)), pch = 16, col = unique(unname(label_vect)), bty = "n")
  # dev.copy(png, file=paste0(folder_name,'/PCA_with_healthy_OG.png'))
  # dev.off()
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
    scale_color_manual(values=levels(pca_data$colors))+
    geom_point(size=4, stroke = 2) +
    xlab(percentVar[1]) +
    ylab(percentVar[2]) +
    theme(legend.position = "none",
          text = element_text(size=20))
  
  ggsave(pca_plot,file=paste0(folder_name,'/PCA_with_healthy.png'))
  
  pca_table<-as.data.frame(arraydata.pca$var.explained*100)
  colnames(pca_table)=c('% variability')
  row.names(pca_table)=paste0('PC', row.names(pca_table))
  write.csv(pca_table,paste0(folder_name,'/PCA_with_healthy_data.csv'))
  
  #Remove healthy patients for PCA without healthy patients
  arraydata.no.healthy<- filtered_array_dat
  for (sample in colnames(filtered_array_dat$E)){
    if (startsWith(sample,'F')==T){
      arraydata.no.healthy <- arraydata.no.healthy[,!colnames(arraydata.no.healthy$E) == sample]
    }
  }
  color_vector<-color_vector[color_vector!='Healthy']
  # pca plot
  pdf(file=paste0(folder_name,'/OG_no_health.pdf'))#Prevent plotting of the MDS, just want the values
  arraydata.no.healthy <- plotMDS(arraydata.no.healthy, pch = 16, gene.selection = "common", col=color_vector)
  dev.off()
  abc <- cmdscale(as.dist(arraydata.no.healthy$distance.matrix), eig=T)
  PC1 <- paste("PC1 (", round(arraydata.no.healthy$var.explained[1]*100, 1), "%)", sep = "")
  PC2 <- paste("PC2 (", round(arraydata.no.healthy$var.explained[2]*100, 1), "%)", sep = "")
  percentVar<-c(PC1,PC2)
  
  # plotMDS(arraydata.no.healthy, pch = 16, gene.selection = "common", col=color_vector,dim.plot=c(1,2))
  # legend(par("usr")[2], par("usr")[4], unique(names(label_vect)), pch = 16, col = unique(unname(label_vect)), bty = "n")
  # dev.copy(png, file=paste0(folder_name,'/PCA_no_healthy_PC1_PC2.png'))
  # dev.off()
  
  label_vect<-label_vect[names(label_vect)!='Healthy']
  
  pca_data<-data.frame(PC1=arraydata.no.healthy$x,PC2=arraydata.no.healthy$y,
                       group=factor(names(label_vect),levels=unique(names(label_vect))),
                       colors=factor(unname(label_vect),levels=unique(unname(label_vect))))
  
  
  
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
    scale_color_manual(values=levels(pca_data$colors))+
    geom_point(size=4, stroke = 2) +
    xlab(percentVar[1]) +
    ylab(percentVar[2]) +
    theme(legend.position = "none",
          text = element_text(size=20))
  
  ggsave(pca_plot,file=paste0(folder_name,'/PCA_no_healthy.png'))
  
  pca_table<-as.data.frame(arraydata.no.healthy$var.explained*100)
  colnames(pca_table)=c('% variability')
  row.names(pca_table)=paste0('PC', row.names(pca_table))
  write.csv(pca_table,paste0(folder_name,'/PCA_no_healthy_data.csv'))
  #
  # plotMDS(arraydata.no.healthy, pch = 16, gene.selection = "common", col=color_vector,dim.plot=c(1,3))
  # legend(par("usr")[2], par("usr")[4], unique(names(label_vect)), pch = 16, col = unique(unname(label_vect)), bty = "n")
  # dev.copy(png, file=paste0(folder_name,'/PCA_no_healthy_PC1_PC3.png'))
  # dev.off()
  #
  # plotMDS(arraydata.no.healthy, pch = 16, gene.selection = "common", col=color_vector,dim.plot=c(2,3))
  # legend(par("usr")[2], par("usr")[4], unique(names(label_vect)), pch = 16, col = unique(unname(label_vect)), bty = "n")
  # dev.copy(png, file=paste0(folder_name,'/PCA_no_healthy_PC2_PC3.png'))
  # dev.off()
  
}


# Plotting functions -------------------



#' Create an MA plot with option of labeling the genes of interest
#'
#' The function creates a 'four quadrant' MA plot, where the FDR and log2FoldChange
#' thresholds dictate the significant up-regulated category, the significant downregulated
#' category, the significant low regulation category and the non-significant category
#' 
#' Genes of interest are labeled in a rectangle for visibility and will have the same 
#' color as the category in which they are in.
#'
#' The plot is created using ggplot2, to save a plot the ggsave() function is recommended
#' It is also recommended to use the following parameters to save the plot.
#' dpi=300 width=21 height=19 units='cm'
#'
#' @param DE_res The differential expression results to be plotted
#' @param genes_of_interest A vector containing gene names to be labelled in the plot
#' To not label any genes, leave as default or provide an empty vector.
#' @param l2FC_thresh The value used as the log2FoldChange threshold
#' @param FDR_thresh The value used to establish significance from adjusted pvalue
#' @param plot_title The title to be give to the plot
#' 
#' @return None
#'
#' @examples
#' 
#' ggsave(paste0("DE_results/",cond_1,"_vs_",cond_2,"/ma_plot.png"),dpi=300,width=21, height=19, units='cm',
#'        plot=maplot_alt(DE_res = DE_gtf, genes_of_interest = genes_of_interest,plot_title = paste0(cond_1," vs ",cond_2))
#' )
#' @export
maplot_alt <- function(DE_res,filter_choice,genes_of_interest=c(),l2FC_thresh=1,FDR_thresh=0.05,plot_title="MA plot"){
  #This ensures that the non-significant appear first on the plot (background) while
  #the significant genes are plotted on the foreground
  DE_res <- DE_res[order(-DE_res[,filter_choice]),]
  DE_res$baseMean_logged <- log2(DE_res$baseMean + 1)
  
  DE_res$Significance <- NA
  sigUp <- which(DE_res$log2FoldChange>l2FC_thresh & DE_res[filter_choice]<FDR_thresh)
  DE_res[ sigUp, "Significance"] <- paste0("up-reg | ",filter_choice,"<",FDR_thresh," (",length(sigUp),")")
  #subset and fill columns based on significant downregulation
  sigDown <- which(DE_res$log2FoldChange<(-l2FC_thresh) & DE_res[filter_choice]<FDR_thresh)
  DE_res[ sigDown, "Significance"] <- paste0("down-reg | ",filter_choice,"<",FDR_thresh," (",length(sigDown),")")
  
  #subset and fill columns based on non-significant upregulation
  nonSig <- which(DE_res[filter_choice]>=FDR_thresh)
  DE_res[ nonSig, "Significance"] <- paste0("non-significant | ",filter_choice,">",FDR_thresh," (",length(nonSig),")")
  
  lowSig <- which(abs(DE_res$log2FoldChange)<=l2FC_thresh & DE_res[filter_choice]<FDR_thresh)
  DE_res[ lowSig, "Significance"] <- paste0("low-difference | ",filter_choice,"<",FDR_thresh," (",length(lowSig),")")
  
  
  labs_data <- DE_res[DE_res$gene_id %in% genes_of_interest, ]
  ggplot(DE_res, aes(x=baseMean_logged, y=log2FoldChange,color=Significance))+
    geom_point(size=0.8)+
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_hline(yintercept=l2FC_thresh, linetype="dashed", color = "black")+
    geom_hline(yintercept=-l2FC_thresh, linetype="dashed", color = "black")+
    xlab(expression('BaseMean')) +
    ylab(expression('log2FoldChange'))+
    ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 1, segment.colour = 'black',show.legend = FALSE,
                              label.size = 0.5) +
    scale_color_manual(breaks = c(paste0("up-reg | ",filter_choice,"<",FDR_thresh," (",length(sigUp),")"),
                                  paste0("down-reg | ",filter_choice,"<",FDR_thresh," (",length(sigDown),")"),
                                  paste0("low-difference | ",filter_choice,"<",FDR_thresh," (",length(lowSig),")"),
                                  paste0("non-significant | ",filter_choice,">",FDR_thresh," (",length(nonSig),")")),
                       values=c("#B31B21","#1465AC","green","darkgray")) +
    guides(color = guide_legend(override.aes = list(size = 7)))+
    
    ggtitle(plot_title) +
    theme_light()+
    theme(text = element_text(size = 10),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 10))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#' Create a volcano plot with option of labeling the genes of interest and/or top
#' significant genes.
#' 
#' The function creates a 'four quadrant' volcano plot, where the FDR and log2FoldChange
#' thresholds dictate the significant up-regulated category, the significant downregulated
#' category, the significant low regulation category and the non-significant category
#' 
#' Genes of interest are labeled in a rectangle for visibility and will have the same 
#' color as the category in which they are in.
#' Top significant genes are also labelled in rectangles but will have black text
#' in order to distinguish them from the genes of interest
#'
#' The plot is created using ggplot2, to save a plot the ggsave() function is recommended
#' It is also recommended to use the following parameters to save the plot.
#' dpi=300 width=21 height=19 units='cm'
#'
#' @param DE_res The differential expression results to be plotted
#' @param genes_of_interest A vector containing gene names to be labelled in the plot
#' To not label any genes, leave as default or provide an empty vector.
#' @param l2FC_thresh The value used as the log2FoldChange threshold
#' @param FDR_thresh The value used to establish significance from adjusted pvalue
#' @param plot_title The title to be give to the plot
#' @param label_top_n The number of genes to label. Genes will be taken in order of
#' significance where the genes with the lowest adjusted p-value are taken first.
#' @param show_non_sig_interestA boolean to indicate if the non-significant genes of
#' interest should be shown.
#' 
#' @return None
#'
#' @examples
#' 
#' ggsave("DE_results/PBMC_2_vs_IgD/volcano_plot_with_top_10.png",dpi=300,width=21, height=19, units='cm',
#'        volcanoplot_alt(DE_res = DE_gtf,genes_of_interest=genes_of_interest,
#'                        plot_title = "PBMC_2_vs_IgD",label_top_n=10)
#' )
#' @export
volcanoplot_alt <- function(DE_res,genes_of_interest=c(),filter_choice,FDR_thresh=0.05, l2FC_thresh=1,plot_title="Volcano plot",label_top_n=0,show_non_sig_interest=TRUE){
  
  
  
  DE_res$Significance <- NA
  sigUp <- which(DE_res$log2FoldChange>l2FC_thresh & DE_res[filter_choice]<FDR_thresh)
  DE_res[ sigUp, "Significance"] <- paste0("up-reg with ",filter_choice,"<",FDR_thresh," (",length(sigUp),")")
  #subset and fill columns based on significant downregulation
  sigDown <- which(DE_res$log2FoldChange<(-l2FC_thresh) & DE_res[filter_choice]<FDR_thresh)
  DE_res[ sigDown, "Significance"] <- paste0("down-reg with ",filter_choice,"<",FDR_thresh," (",length(sigDown),")")
  
  #subset and fill columns based on non-significant upregulation
  nonSig <- which(DE_res[filter_choice]>=FDR_thresh)
  DE_res[ nonSig, "Significance"] <- paste0("non-significant with ",filter_choice,">",FDR_thresh," (",length(nonSig),")")
  
  lowSig <- which(abs(DE_res$log2FoldChange)<=l2FC_thresh & DE_res[filter_choice]<FDR_thresh)
  DE_res[ lowSig, "Significance"] <- paste0("low-regulation with ",filter_choice,"<",FDR_thresh," (",length(lowSig),")")
  
  #Find the lowest significant point to draw a horizontal ax-line at that point on the axis
  #Sort using the filter choice
  DE_res <- DE_res[order(DE_res[[filter_choice]]),]
  #Find the 'first' non significant gene
  gene_for_sig_line<-DE_res$gene_id[DE_res$Significance==paste0("non-significant with ",filter_choice,">",FDR_thresh," (",length(nonSig),")")][1]
  #Create a dataframe of necessary columns and convert pvalue to -log10 (as is it's form in the plot)
  find_lowest_sig<-DE_res[,c('gene_id','pvalue')]
  find_lowest_sig$pvalue<- -log10(find_lowest_sig$pvalue)
  #Extract value based on the gene.
  lowest_sig<-find_lowest_sig$pvalue[find_lowest_sig$gene_id==gene_for_sig_line]
  labs_data <- DE_res[DE_res$gene_id %in% genes_of_interest, ]
  if (show_non_sig_interest==FALSE){
    labs_data <- labs_data[labs_data[filter_choice]<FDR_thresh & abs(labs_data$log2FoldChange)>l2FC_thresh,]
  }
  if (label_top_n>0){
    top_res <- DE_res
    labs_data_top <- top_res[1:label_top_n,]
  }else{
    labs_data_top <- DE_res[FALSE,]
  }
  ggplot(DE_res, aes(x = log2FoldChange, y = -log10(pvalue), color=Significance)) +
    geom_point(size=0.8) +
    geom_hline(yintercept=lowest_sig, linetype="dashed", color = "black")+
    xlab(expression('Log'[2]*'Fold Change')) +
    ylab(expression('-log10(pvalue)')) +
    geom_vline(xintercept = l2FC_thresh, linetype="dashed", color = "black") +
    geom_vline(xintercept = -l2FC_thresh, linetype="dashed", color = "black") +
    ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 1, segment.colour = 'black',show.legend = FALSE,
                              label.size = 1)+
    
    ggrepel::geom_label_repel(data = labs_data_top, mapping = aes(label = gene_id),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 0.5, colour = 'black',show.legend = FALSE) +
    
    scale_color_manual(breaks = c(paste0("up-reg with ",filter_choice,"<",FDR_thresh," (",length(sigUp),")"),
                                  paste0("down-reg with ",filter_choice,"<",FDR_thresh," (",length(sigDown),")"),
                                  paste0("non-significant with ",filter_choice,">",FDR_thresh," (",length(nonSig),")"),
                                  paste0("low-regulation with ",filter_choice,"<",FDR_thresh," (",length(lowSig),")")),
                       values=c("#B31B21","#1465AC","darkgray","green")) +
    guides(color = guide_legend(override.aes = list(size = 7)))+
    ggtitle(plot_title) +
    theme_light()+
    theme(text = element_text(size = 10),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 10))
  
}


# Unused functions -------------------
#' Annotation function
#' 
#' This function updates the gene segment of the Elist object using
#' annotation data from agilent for the specific microarray used.
#' 
#' @param data The E list object to be updated
#' @param annotation_file A dataframe of the anontation file
#' 
#' @return data The annotated Elist object
#'
#' @examples
#'
#' @export
perform_annotation <- function(data,annotation_file){
  annot_cols_interest <- c('ProbeID','GeneSymbol')
  annotation_file <- annotation_file[,annot_cols_interest]
  colnames(annotation_file)=c('ProbeName','gene_id')
  
  
  my_gene_col<-data[,c('ProbeName','SystematicName')]
  
  my_gene_col <- merge(my_gene_col,annotation_file,by='ProbeName')
  my_gene_col <- my_gene_col[match(data$ProbeName, my_gene_col$ProbeName),]
  
  data$gene_id <- my_gene_col$gene_id
  
  return(data)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


create_heat_genes_interest <- function(exp_list,genes_of_interest,name_file){
  saveWidget(
    plot_l2fc_heatmap(experiment_list = lst, genes_of_interest = genes_of_interest,cluster=T,name_file ),
    file="DE_results/heatmap_L2FC_summary_clustered.html")
  #Savewidget creates a library folder (for some reason) so we remove it with unlist
  unlink(paste0("DE_results/heatmap_L2FC_summary_clustered","_files"), recursive=TRUE)
  
  saveWidget(
    plot_l2fc_heatmap(experiment_list = lst, genes_of_interest = genes_of_interest,cluster=F,name_file ),
    file="DE_results/heatmap_L2FC_summary_no_cluster.html")
  #Savewidget creates a library folder (for some reason) so we remove it with unlist
  unlink(paste0("DE_results/heatmap_L2FC_summary_no_cluster","_files"), recursive=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' A hardcoded function to create subsets and design a matrix from those subsets
#' 
#' The function is hardcoded for a specific file with specific headers. It creates
#' pre-determined subsets within this file.
#' Each subset is saved as a separate csv file to enable confirmation of it's validity
#' by a clinician
#' 
#' 
#' 
#' @param patient_sample_dta The file which will be subsetted
#' @param save_path The path where the subsets are saved
#' 
#' @return my_mm the resulting model matrix using the subsets
#'
#' @examples
#'
#' @export
create_custom_CRC_matrix <- function(patient_sample_dta,save_path){
  #Creating my own matrix
  ####
  sampdta_for_matrix <- patient_sample_dta[,c('SampleID','progressor_status','inflammation_status','neoplasia_status','kryptDist_status')]
  samp_dta_for_print <- sampdta_for_matrix
  #Establish the subsets
  NP_inflammation <- samp_dta_for_print[samp_dta_for_print$progressor_status=='NP' & 
                                          samp_dta_for_print$inflammation_status!=0,]
  sampdta_for_matrix$NP_inflammation <- 0
  sampdta_for_matrix$NP_inflammation[sampdta_for_matrix$SampleID %in% NP_inflammation$SampleID] <- 1
  
  write.csv(NP_inflammation,file = paste0(save_path,'/NP_inflammation.csv'),row.names = F)
  
  ####
  
  P_inflammation_non_neo <- samp_dta_for_print[samp_dta_for_print$progressor_status=='P' & 
                                                 samp_dta_for_print$inflammation_status!=0 &
                                                 samp_dta_for_print$neoplasia_status==0,]
  sampdta_for_matrix$P_inflammation_non_neo <- 0
  sampdta_for_matrix$P_inflammation_non_neo[sampdta_for_matrix$SampleID %in% P_inflammation_non_neo$SampleID] <- 1
  
  write.csv(P_inflammation_non_neo,file = paste0(save_path,'/P_inflammation_non_neo.csv'),row.names = F)
  
  ####
  
  NP_non_inflammation <- samp_dta_for_print[samp_dta_for_print$progressor_status=='NP' & 
                                              samp_dta_for_print$inflammation_status==0,]
  sampdta_for_matrix$NP_non_inflammation <- 0
  sampdta_for_matrix$NP_non_inflammation[sampdta_for_matrix$SampleID %in% NP_non_inflammation$SampleID] <- 1
  
  write.csv(NP_non_inflammation,file = paste0(save_path,'/NP_non_inflammation.csv'),row.names = F)
  
  ####
  
  P_non_inflammation_non_neo <- samp_dta_for_print[samp_dta_for_print$progressor_status=='P' & 
                                                     samp_dta_for_print$inflammation_status==0 & 
                                                     samp_dta_for_print$neoplasia_status==0,]
  sampdta_for_matrix$P_non_inflammation_non_neo <- 0
  sampdta_for_matrix$P_non_inflammation_non_neo[sampdta_for_matrix$SampleID %in% P_non_inflammation_non_neo$SampleID] <- 1
  
  write.csv(P_non_inflammation_non_neo,file = paste0(save_path,'/P_non_inflammation_non_neo.csv'),row.names = F)
  
  ####
  
  P_neoplasia <-samp_dta_for_print[samp_dta_for_print$progressor_status=='P' & 
                                     samp_dta_for_print$neoplasia_status>=2 & 
                                     samp_dta_for_print$neoplasia_status!=6,]
  sampdta_for_matrix$P_neoplasia <- 0
  sampdta_for_matrix$P_neoplasia[sampdta_for_matrix$SampleID %in% P_neoplasia$SampleID] <- 1
  
  
  write.csv(P_neoplasia,file = paste0(save_path,'/P_neoplasia.csv'),row.names = F)
  
  ####
  
  NP_normal <- samp_dta_for_print[samp_dta_for_print$progressor_status=='NP' &
                                    samp_dta_for_print$inflammation_status==0 &
                                    samp_dta_for_print$neoplasia_status==0 &
                                    samp_dta_for_print$kryptDist_status==0,]
  sampdta_for_matrix$NP_normal <- 0
  sampdta_for_matrix$NP_normal[sampdta_for_matrix$SampleID %in% NP_normal$SampleID] <- 1
  
  write.csv(NP_normal,file = paste0(save_path,'/NP_normal.csv'),row.names = F)
  
  ####
  
  P_normal <- samp_dta_for_print[samp_dta_for_print$progressor_status=='P' &
                                   samp_dta_for_print$inflammation_status==0 &
                                   samp_dta_for_print$neoplasia_status==0 &
                                   samp_dta_for_print$kryptDist_status==0,]
  sampdta_for_matrix$P_normal <- 0
  sampdta_for_matrix$P_normal[sampdta_for_matrix$SampleID %in% P_normal$SampleID] <- 1
  
  write.csv(P_normal,file = paste0(save_path,'/P_normal.csv'),row.names = F)
  
  
  sampdta_for_matrix <- sampdta_for_matrix[,6:12]
  sampdta_for_matrix <- sampdta_for_matrix[,c('P_inflammation_non_neo','NP_normal')]
  
  # my_mm <- model.matrix(~ NP_inflammation + P_inflammation_non_neo + NP_non_inflammation +
  #                         P_non_inflammation_non_neo + P_neoplasia + NP_normal + P_normal,sampdta_for_matrix)
  
  my_mm <- model.matrix(~ P_inflammation_non_neo + NP_normal,sampdta_for_matrix)
  
  colnames(my_mm) <- make.names(colnames(my_mm))
  
  return(my_mm)
}



#' A hardcoded function to create contrasts for differential gene expression
#' 
#' The function is hardcoded to produce the comparisons that we want to perform
#' with the subsets created from the create_custom_CRC_matrix function
#' 
#' 
#' 
#' @param the_mm The model matrix to be used as levels
#' 
#' @return my_mm the resulting model matrix using the subsets
#'
#' @examples
#'
#' @export
prepare_contrasts_for_CRC <- function(the_mm){
  ctrst <- makeContrasts(#NP_inflammation-P_inflammation_non_neo,
    # NP_non_inflammation-P_non_inflammation_non_neo,
    # NP_non_inflammation-NP_inflammation,
    # P_non_inflammation_non_neo-P_inflammation_non_neo,
    # NP_non_inflammation-P_inflammation_non_neo,
    # NP_inflammation-P_non_inflammation_non_neo,
    # P_non_inflammation_non_neo-P_neoplasia,
    # P_inflammation_non_neo-P_neoplasia,
    NP_normal-P_inflammation_non_neo,
    levels=the_mm)
  colnames(ctrst) <- c(#'NP.inflammation_vs_P.inflammation.non.neo',
    # 'NP.non.inflammation_vs_P.non.inflammation.non.neo',
    # 'NP.non.inflammation_vs_NP.inflammation',
    # 'P.non.inflammation.non.neo_vs_P.inflammation.non.neo',
    # 'NP.non.inflammation-P.inflammation.non.neo',
    # 'NP.inflammation_vs_P.non.inflammation.non.neo',
    # 'P.non.inflammation.non.neo_vs_P.neoplasia',
    # 'P.inflammation.non.neo_vs_P.neoplasia',
    'NP.normal_vs_P.inflammation.non.neo')
  
  
  
  # ctrst <- makeContrasts(NP_non_inflammation-NP_inflammation,
  #                        P_norm_jonas-P_inflammation_non_neo,
  #                            levels=the_mm)
  # colnames(ctrst) <- c('NP.non.inflammation_vs_NP_non_inflammation',
  #                      'P.norm.jonas_vs_P.inflammation')
  
  return(ctrst)
}
