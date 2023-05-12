library(WGCNA)
library(ggrepel)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)





# pipeline functions  -------------------

#' Function which obtains and formats the data needed for the WGCNA analysis
#'
#' The function will automatically search for the required csv files
#' based on the inputted experiment. It then formats the data in accordance
#' to the requirements of the WGCNA package
#' 
#' @param targeted_experiment The name of the experiment to be read and formatted
#' @param genes_of_interest Adds the genes of interest as the clinical data 
#' 
#' @return return_list A list containing the raw data and the clinical data
#'
#' @examples
#' target_exp <- IgM.acD40.IL4.IL21_vs_PBMC.1"
#' genes_of_interest <- c("AICDA","APOBEC1","APOBEC2","APOBEC3A","APOBEC3B",
#' "APOBEC3C","APOBEC3D","APOBEC3E","APOBEC3F","APOBEC3G",
#' "APOBEC3H","APOBEC4")
#' returned_list <- prep_WGCNA_data(target_exp,genes_of_interest)
#'
#' @export
prep_WGCNA_data <- function(targeted_experiment, genes_of_interest){
  #Declare that we want to use this specific differential expression file
  dir.create(paste0("WGCNA_results/",target_exp))
  
  
  de_raw<-read.csv(paste0("DE_results/",target_exp,"/DE_raw_data.csv"))
  de_raw = subset(de_raw, select = -c(2:7) )
  
  # de_gtf<-read.csv(paste0("DE_results/",target_exp,"/DE_results_with_gtf.csv"))
  #Prep our clinical data
  
  if (is.null(genes_of_interest)==F){
    datTraits<-de_raw[de_raw$gene_id %in% genes_of_interest,]
    rownames(datTraits)=datTraits[,1]
    datTraits<-datTraits[,-1]
    datTraits_t<-transpose(datTraits)
    rownames(datTraits_t)=colnames(datTraits)
    colnames(datTraits_t)=rownames(datTraits)
  }else{
    datTraits <- de_raw[FALSE,]
    datTraits_t<-datTraits[,-1]
  }
  #Subset for testing purposes
  # de_gtf<-de_gtf[1:3500,]
  
  gene_list<-de_raw$gene_id
  
  de_raw<-de_raw[de_raw$gene_id %in% gene_list,]
  
  rownames(de_raw)=de_raw$gene_id
  de_raw<-de_raw[,-1]
  
  de_raw_t<-transpose(de_raw)
  rownames(de_raw_t)=colnames(de_raw)
  colnames(de_raw_t)=rownames(de_raw)
  de_raw<-de_raw_t
  
  return_list <- list(de_raw,datTraits_t)
  return(return_list)
}


prep_WGCNA_data_V2<-function(count_folder,datTraits){
  #Check if samples are present
  count_files<-list.files(count_folder)
  count_files<-unlist(strsplit(count_files,'.counts'))
  missing_samples<-datTraits$sample[!datTraits$sample %in% count_files]
  if (length(missing_samples)>0){
    message('The following samples have not been found and will be removed')
    message(missing_samples)
    
    datTraits<-datTraits[!datTraits$sample %in% missing_samples,]
  }
  
  main_matrix<-data.frame(NULL)
  for(counts in datTraits$sample){
    temp_df<-read.delim(paste0(count_folder,'/',counts,'.counts'),header=F,row.names=1)
    temp_df_t<-t(temp_df)
    row.names(temp_df_t)=counts
    colnames(temp_df_t)=row.names(temp_df)
    
    if(nrow(main_matrix)==0){
      main_matrix<-temp_df_t
    }else{
      main_matrix<-rbind(main_matrix,temp_df_t)
    }
    
  }
  
  
  #Convert all columns of datTraits to numberic (except sample)
  cols_to_keep<-as.vector(colnames(datTraits)[colnames(datTraits)!='sample'])
  datTraits_sub<-datTraits[,cols_to_keep,drop=F]
  row.names(datTraits_sub)=datTraits$sample
  
  for(col in colnames(datTraits_sub)){
    datTraits_sub[[col]]<-as.numeric(datTraits_sub[[col]])
  }

  return(list(main_matrix,datTraits_sub))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which preforms some quality control on the dataset
#'
#' The function checks for outlier samples and removes samples with too few values
#' An indication/print is given in the event that samples had to be removed.
#' 
#' @param datExpr0 The WGCNA formatted data.
#' 
#' @return datExpr0 The WGCNA data with the removed sample if applicable
#'
#' @examples
#' datExpr0 <- pre_WGCNA_QC(datExpr0)
#'
#' @export
pre_WGCNA_QC <- function(datExpr0){
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  return(datExpr0)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which determines the cluster 'under the line'
#'
#' The function serves as a means to remove samples which are outliers on the dendrograme
#' By inputting a cutoff height, the function will remove that sample and re-cluster
#' the remaining ones. The function also selectes the cluster to keep, by default it
#' is 0. This step should be double checked by users before proceeding to the next script.
#' 
#' If a cutoff is implemented, the cluster required may no longer be 0. The clusters
#' can be examined using the table(clust) command, which is commented out below the
#' if statement
#' 
#' @param sample_tree The previously clustered tree
#' @param datExpr0 The WGCNA data
#' @param cutoff_inputted The cutoff height. By default it is null so no cutoff
#' @param cluster_chosen The cluster to be kept and used downstream.
#' 
#' @return return_list A list containing the updated sampleTree, the colors associated
#' to the tree and the subsetted WGCNA data to only contain the required samples.
#'
#' @examples
#' sampleTree = hclust(dist(datExpr0), method = "average")
#' returned_list <- cluster_under_line(sample_tree,datExpr0)
#'
#' @export
cluster_under_line <- function(sample_tree,datExpr0,datTraits,cutoff_inputted=NULL,cluster_chosen=0){
  
  if (is.null(cutoff_inputted)==F){
    clust = cutreeStatic(sampleTree, cutHeight = cutoff_inputted, minSize = 10)
  }
  clust = cutreeStatic(sampleTree, minSize = 10)
  # table(clust)
  
  keepSamples = (clust==cluster_chosen)
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE)
  
  return_list<-list(sampleTree2,traitColors,datExpr)
  return(return_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Creates a list of genes for visEAGO using WGCNA results
#'
#' Given a list of modules/colors of interest, the function will open the first 
#' 'gene_info' file it finds in the WGCNA folder. Using the modules of interest
#' provided it will create a list containing the genes within each module/color
#' Each element of the list is named after it's respective module.
#' A background element is also created, it contains all the ENTREZID in the genome
#' genome according to the annotationdbi package.
#' 
#' The function converts the gene symbols found into ENTREZID
#' 
#'
#' @param modules_of_interest a vector of modules of interest
#' @param organism the annotation dbi object used to convert gene symbols to ENTREZID
#' @param target_exp the name of the folder with the results to process
#' @param convert Boolean to check if conversion between gene symbol and entrez id is necessary
#'
#' @return my_gene_list A named list containing gene names for each submitted differential expression
#'
#' @examples
#'
#' @export
prep_WGCNA_for_overrepresentation <-function(modules_of_interest,target_exp,organism,convert=T){
  #Find and extract gene info results from WGCNA results
  list_of_files<-list.files(paste0('WGCNA_results/',target_exp,'/coexpression_genes'))
  target_files <- list_of_files[grep('_geneInfo.csv', list_of_files)[1]]
  gene_info_path<-paste0('WGCNA_results/',target_exp,'/coexpression_genes/',target_files)
  gene_info<-read.csv(gene_info_path)
  
  if(modules_of_interest=='All'){
    modules_of_interest<-unique(gene_info$moduleColor)
  }
  
  
  #Prep visEAGO gene list based on WGCNA module results
  my_gene_list <-list()
  if (convert==T){
    my_gene_list[['background']]<-keys(organism,keytype = 'ENTREZID')
  }else{
    my_gene_list[['background']]<-keys(organism,keytype = 'SYMBOL')
  }
  
  for(color in unique(gene_info$moduleColor)){
    if(color %in% modules_of_interest){
      my_gene_list[[color]]<-gene_info$geneSymbol[gene_info$moduleColor==color]
      if(convert==T){
        my_gene_list[[color]]<-as.vector(mapIds(organism, keys = my_gene_list[[color]], keytype = "SYMBOL", column="ENTREZID"))
        my_gene_list[[color]] <- my_gene_list[[color]][!is.na(my_gene_list[[color]])]
      }
    }
  }
  return(my_gene_list)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which prepares a correlation matrix using spearmans correlation
#' along with the two.sided alternative hypothesis.
#' 
#' The function also creates a formatted pvalue matrix for a diagonal correlation
#' plot. 
#' Both the correlation matrix and the diagonal pvaue matrix are returned
#' 
#'
#' @param mat The dataframe containing the clinical data
#' @param ... Any parameter to be given to cor.test function
#'
#' @return return_list A list containing the diagonal formatted pvalue matrix and the
#' correlation matrix
#'
#' @examples
#'
#' @export
calculate_corr_and_pval_matrix <- function(mat,...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  corr.mat<- matrix(NA, n, n)
  diag(corr.mat) <- 0
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method='spearman',alternative='t',exact = FALSE,...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      corr.mat[i, j] <- corr.mat[j, i] <- tmp$estimate
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat_for_format<-p.mat
  for (idx in 1:nrow(p.mat_for_format)){
    if (idx==1){
      p.mat_for_format[idx,idx]<-''
    }else{
      p.mat_for_format[idx,1:idx]<-''
    }
  }
  colnames(corr.mat)<-colnames(p.mat)
  row.names(corr.mat)<-row.names(p.mat)
  
  return_list<-list(p.mat_for_format,corr.mat)
  return(return_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' The function converts the pvalue matrix from \code{calculate_corr_and_pval_matrix} to
#' a matrix with the scientific notation equivalent of the pvalue matrix.
#' The function then reverses the order enabling the use of the formatted matrix in
#' the \code{plot_clinical_corr_plot} for the plotting of the correlation matrix.
#' Any non-significant pvalues are replaced with blank spaces.
#' 
#'
#' @param p.mat The diagonal pvalue matrix for the correlation values
#'
#' @return formatted_p The formatted pvalue matrix
#'
#' @examples
#'
#' @export
format_p_val_matrix <-function(p.mat){
  formatted_p<-signif(as.numeric(p.mat), digits=3)

  for (idx in 1:length(formatted_p)){
    if (is.na(formatted_p[idx])==F){
      if (formatted_p[idx]>0.05){
        formatted_p[idx]<-NA
      }
    }
  }
  formatted_p<-formatC(formatted_p, format = "e", digits = 2)
  for (idx in 1:length(formatted_p)){
    if (formatted_p[idx]==" NA"){
      formatted_p[idx]<-''
    }
  }
  
  formatted_p<-rev(formatted_p)
  
  return(formatted_p)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Creates a dataframe of GO's within modules of interest, can save gprofiler results
#'
#' The function takes in a list of modules of interest, a result location, ontology,
#' and organism. It then performs a Gprofiler over-representation for each module (individually)
#' and stores the GO results in a dataframe. The dataframe is returned
#' 
#' If save_path is not null, the gprofiler results will be saved to the designated
#' location in both csv and html (interactive plot) format.
#' 
#'
#' @param modules_of_interest a vector of modules of interest
#' @param path_to_modules The folder path to the gene list csv files for the module 
#' @param ontology the ontology that will be returned in the dataframe (ex: GO:BP)
#' @param organism organism to use with gprofiler analysis (ex: hsapien)
#' @param save_path The folder path to save results if gprofiler results are to be saved
#'
#' @return module_GO A dataframe containing all the GOs (ID) found, their module, and the term name
#'
#' @examples
#'
#' @export
gprofiler_modules_interest<-function(modules_of_interest,path_to_modules,ontology,organism,save_path=NULL){
  #If a save_path was given, create folders to save the data
  if (is.null(save_path)==FALSE){
    dir.create(paste0(save_path,'/gprofiler_results'))
    save_fig_location<-paste0(save_path,'/gprofiler_results/figures')
    save_data_location<-paste0(save_path,'/gprofiler_results/data_files')
    dir.create(save_fig_location)
    dir.create(save_data_location)
  }
  
  module_GO<-as.data.frame(NULL)
  for (module in modules_of_interest){
    #Read and format the genes from the module of interest
    gene_vect<-read.csv(paste0(path_to_modules,'/',module,'.csv'))
    gene_vect<-as.vector(gene_vect$geneSymbol)

    gostres <- gost(query = gene_vect,organism = organism)
    #Checks if the query returned any results
    if (is.null(gostres)==F){
      #store figure results
      p <- gostplot(gostres, capped = T, interactive = T)
      
      #Produce data results - conversion of parent column (collapse list)
      gprof_dta<-gostres$result
      for (idx in 1:nrow(gprof_dta)){
        gprof_dta$parents[idx]<-paste(unlist(gprof_dta$parents[idx]), collapse = '/')
      }
      gprof_dta$parents<-unlist(gprof_dta$parents)
      
      #Add to list for subsequent plotting
      if (ontology %in% gprof_dta$source == TRUE){
        extracted_results<-gprof_dta[c('term_id','term_name','p_value')][gprof_dta$source==ontology,]
        extracted_results$module_col<-rep(strsplit(module,'_')[[1]][1],nrow(extracted_results))
        module_GO<-rbind(module_GO,extracted_results)
      }
      
      #If a save_path was given, save the data
      if (is.null(save_path)==FALSE){
        module_name<-strsplit(module,'.csv')[[1]]
        save_name_fig<-paste0(save_fig_location,'/',module_name,'_overview.html')
        save_name_data<-paste0(save_data_location,'/',module_name,'_data.csv')
        #save plot
        saveWidget(p,file=save_name_fig)
        #Save data
        write.csv(x=gprof_dta,file=save_name_data,quote=FALSE)
        #Delete extra folder created by htmlwidgets
        unlink(paste0(save_fig_location,'/',module_name,'_overview_files'),recursive = T)
      }
    }
  }
  module_GO<-as.data.frame(module_GO)
  return (module_GO)
}



# plotting functions  -------------------


#' Function which plots the possible power selection for an experiment
#'
#' The function takes in WGCNA data and calculates the soft thresholding power. 
#' The user can then select the proper power for  a specific experiment based on 
#' these plots. The plot also shows a red line at the 0.8 height, users should 
#' select the power value above that line where the values start plateau-ing.
#' 
#' @param datExpr The WGCNA formatted data
#' 
#' @examples
#' sampleTree = hclust(dist(datExpr0), method = "average")
#' returned_list <- cluster_under_line(sample_tree,datExpr0)
#'
#' @export
plot_power <-function(datExpr){
  #power settings
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr,power=powers,verbose=5)#, power=powers, verbose = 5)
  
  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.80,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which calculates and plots the distribution of the modules in 
#' accordance to the clinical data
#'
#' The function starts by qantifying each module trait association. It determines
#' which modules are significanly associated with clinical traits. It does so
#' by correlating eigengene with extranal traits and looks for the significant associations
#' 
#' Based on these associations, a pvalue can be given, and therefore the significance
#' of each module in regards to the clinical value can be illustrated on the heatmap.
#' 
#' @param net The object resulting from the \code{blockwiseModules()}
#' @param datExpr The WGCNA formatted data
#' @param datTraits The WGCNA formatted clinical data
#' 
#' @examples
#' load("WGCNA_results/IgM.acD40.IL4.IL21_vs_PBMC.1/IgM.acD40.IL4.IL21_vs_PBMC.1_part1.RData")
#' 
#' net = blockwiseModules(datExpr, power = target_power,
#'                        TOMType = "unsigned", minModuleSize = 30,
#'                        reassignThreshold = 0, mergeCutHeight = 0.25,
#'                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#'                        saveTOMs = F,
#'                        maxBlockSize = 20000,
#'                        verbose = 3,nThreads=7)
#'                        
#' plot_labeled_heatmap(net,datExpr,datTraits)
#' 
#' @export
plot_labeled_heatmap <- function(moduleTraitCor,moduleTraitPvalue,datTraits,MEs){
  # # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colorMatrix = NULL,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which creates the csv results for a single clinical parameter
#' 
#' The function, based on the inputted clinical data will create two csv files
#' One which contains the full gene/statistical results in association with that
#' clinical value.
#' Another which is the gene list of the genes deemed significant based on their 
#' pvalue (<0.05)
#' 
#' @param clin_dat The clinical data type to use
#' @param datTraits The WGCNA formatted clinical data
#' @param datExpr The WGCNA formatted data
#' @param file_location
#' 
#' @examples
#' 
#' csv_res_single_clinical("AICDA",block_modules=net,datTraits,datExpr)
#' 
#' @export
csv_res_single_clinical <- function(clin_dat,datTraits,datExpr,MEs,nSamples,moduleColors,gene_names=NULL,file_location='WGCNA_results'){
  
  if (clin_dat %in% colnames(datTraits)==F){
    print("Clinical type selected not in data")
    return()
  }
  # Define variable main_clinical containing the clinical column of interest
  main_clinical = as.data.frame(datTraits[clin_dat])
  names(main_clinical) = clin_dat
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(cor(datExpr, main_clinical, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(main_clinical), sep="")
  names(GSPvalue) = paste("p.GS.", names(main_clinical), sep="")
  if(is.null(gene_names)==F){
    geneInfo0 = data.frame(geneSymbol = gene_names,
                           moduleColor = moduleColors,
                           geneTraitSignificance,
                           GSPvalue)
  }else{
    geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
                           moduleColor = moduleColors,
                           geneTraitSignificance,
                           GSPvalue)
  }
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, main_clinical, use = "p")));
  
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  clin_dat<-gsub(x = clin_dat,pattern = '-',replacement = '\\.')
  GS_clin_dat <- paste0('GS.',clin_dat)
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[GS_clin_dat]))
  
  geneInfo = geneInfo0[geneOrder, ]
  p_GS_clin_dat <- paste0('p.',GS_clin_dat)
  sig_genes <- geneInfo[geneInfo[p_GS_clin_dat]<0.05,]
  sig_genes <- sig_genes[order(sig_genes[p_GS_clin_dat]),]
  

  write.csv(geneInfo, file = paste0(file_location,"/coexpression_genes/",clin_dat,"_geneInfo.csv"))
  write.csv(sig_genes, file = paste0(file_location,"/coexpression_genes/",clin_dat,"_sig_genes.csv"))


  
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which plots a standard heatmap showing the intensity of the clinical
#' data for each patient within the WGCNA analysis
#' 
#' 
#' @param my_traits The clinical data in dataframe format
#' @param my_colors A numbers2colors object, as returned by \code{cluster_under_line}
#' @param sort_by The WGCNA formatted data
#' 
#' @examples
#' 
#' @export
plot_clinical_heatmap<-function(my_traits,my_colors,sort_by){
  
  my_table<-as.data.frame(t(my_colors[order(-my_traits[[sort_by]]),]))
  colnames(my_table)=row.names(my_traits[order(-my_traits[[sort_by]]),])
  
  my_table %>%
    tibble::rownames_to_column() %>%
    mutate(rowname = as.numeric(rowname)) %>%
    pivot_longer(-1) %>%
    mutate(name = factor(name, levels = unique(name)))%>%
    ggplot(aes(x = name, y = rowname, fill = value)) +
    xlab("Patient IDs") + 
    ylab("Clinical values")+
    geom_tile(color = "gray50") +
    scale_fill_identity() +
    scale_x_discrete(position = "bottom",
                     expand = c(0, 0)) +
    scale_y_reverse(breaks = seq(nrow(my_table)),
                    labels=colnames(my_traits),
                    expand = c(0, 0))+
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90,vjust=0.5,hjust=-5,size=12),
      axis.text.y = element_text(size=12))+
    coord_equal()
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Function which plots a correlation matrix for the clinical data
#' 
#' The plot is automatically saved to 'WGCNA_results' folder
#' 
#' @param my_traits Dataframe of clinical values
#' @param correlation_matrix The correlation matrix as calculated from \code{calculate_corr_and_pval_matrix}
#' @param formatted_pvals a pvalue matrix as created by \code{format_p_val_matrix}
#'
#' @examples
#'
#' @export
plot_clinical_corr_plot<-function(my_traits,correlation_matrix,formatted_pvals,file_location='WGCNA_results'){
  col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))
  #Add the text
  pos <- expand.grid(1:ncol(my_traits), ncol(my_traits):1)
  #Plot the correlation plot
  png(paste0(file_location,"/clin_correlation_matrix.png"),width = 1250,height = 1250)
  corrplot(correlation_matrix, method="color", col=col(100),
                 type="upper",
                 # addCoef.col = "black", # Add coefficient of correlation
                 tl.col="black", tl.srt=45, #Text label color and rotation
                 # hide correlation coefficient on the principal diagonal
                 diag=FALSE
  )
  text(rev(pos), formatted_pvals)
  dev.off()
}


#' Function which plots a correlation matrix for the clinical data
#' 
#' The plot is automatically saved to 'WGCNA_results' folder
#' 
#' @param my_traits Dataframe of clinical values
#' @param correlation_matrix The correlation matrix as calculated from \code{calculate_corr_and_pval_matrix}
#' @param formatted_pvals a pvalue matrix as created by \code{format_p_val_matrix}
#'
#' @examples
#'
#' @export
plot_clinical_corr_plot_V2<-function(my_traits,correlation_matrix,formatted_pvals,file_location='WGCNA_results'){
  col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))
  #Add the text
  pos <- expand.grid(1:ncol(my_traits), ncol(my_traits):1)
  #Plot the correlation plot
  svg(paste0(file_location,"/clin_correlation_matrix.svg"),width = 10,height = 10)
  corrplot(correlation_matrix, method="color", col=col(100),
           type="upper",
           sig.level = 0.05,
           # addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           # hide correlation coefficient on the principal diagonal
           diag=FALSE
  )
  text(rev(pos), formatted_pvals)
  dev.off()
}


