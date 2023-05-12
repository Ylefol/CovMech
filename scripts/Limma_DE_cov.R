source('R/my_limma_functions.R')

load('data/processed_microarray_data.rdata')

#Re-order the rows -- is this really necessary? Maybe not, but it can't hurt either
patient_sample_dta <- patient_sample_dta[match(colnames(rna_biop_dataa$E), patient_sample_dta$SampleID),]
dir.create('DE_results')
dir.create('DE_results/Robjects_checkpoints')
save(rna_biop_dataa,patient_sample_dta, file='DE_results/Robjects_checkpoints/final_data_with_patient_data.EList.rdata', version = 2)

lst <- c('critical.cov.D1_vs_non_critical.cov.D1',
         'all.cov.D1_vs_healthy.controls',
         'critical.cov.D3_vs_non_critical.cov.D3',
         'all.cov.D3_vs_healthy.controls',
         'critical.cov.D8_vs_non_critical.cov.D8',
         'all.cov.D8_vs_healthy.controls',
         'critical_cov_D1_vs_healthy.controls',
         'non_critical_cov_D1_vs_healthy.controls',
         'critical_cov_D3_vs_healthy.controls',
         'non_critical_cov_D3_vs_healthy.controls',
         'critical_cov_D8_vs_healthy.controls',
         'non_critical_cov_D8_vs_healthy.controls')

#Viremia comparisons
lst <- c('viremia.cov.D1_vs_non_viremia.cov.D1',
         'viremia.cov.D3_vs_non_viremia.cov.D3',
         'viremia.cov.D8_vs_non_viremia.cov.D8',
         'viremia_cov_D1_vs_healthy.controls',
         'non_viremia_cov_D1_vs_healthy.controls',
         'viremia_cov_D3_vs_healthy.controls',
         'non_viremia_cov_D3_vs_healthy.controls',
         'viremia_cov_D8_vs_healthy.controls',
         'non_viremia_cov_D8_vs_healthy.controls')


#Location of subsets to be used
subset_location <- 'data/subsets_to_use/'
targets<-readTargets("data/micro_arr_metadata.txt", row.names="ID")

#Subset target file and biop data
dir.create('my_PCAs')
target_param<-'Severe'
target_param<-'Viremi'
for (comp in lst){
  return_list<-filter_data_based_on_subsets(rna_biop_dataa,targets,comp,subset_location)
  new_targets<-return_list[[1]]
  new_biop_dta<-return_list[[2]]
  #Set colors
  new_targets[[target_param]][new_targets[[target_param]] == ""] <- "Healthy"
  colors.param <- rep("#4daf4a", nrow(new_targets))
  colors.param[new_targets[[target_param]] == "Y"] <- "#e41a1c"
  colors.param[new_targets[[target_param]] == "N"] <- "#377eb8"
  
  #Set labels
  for_label_param<-colors.param
  names(for_label_param)=new_targets[[target_param]]
  #Run PCA function
  for_save<-paste0('my_PCAs/',comp)
  my_PCA_plotting_function_custom(new_biop_dta,colors.param,for_label_param,for_save)
}


gtf_path <- 'data_files/gencode.v35.annotation.gtf'
gtf_path <- NULL
gene_annotation<-NULL

limma_DE_wrapper(rna_biop_dataa,comparison_list=lst,subset_path=subset_location,sample_dta=patient_sample_dta,
                 gene_annot=gene_annotation,gtf_path=gtf_path,filter_choice='pvalue',genes_interest=c(),species='hg38')




