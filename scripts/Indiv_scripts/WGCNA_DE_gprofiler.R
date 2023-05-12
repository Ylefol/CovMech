

library(gprofiler2)
library(htmlwidgets)
setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")

my_target_dir<-'WGCNA_results/WGCNA_module_gene_lists'
found_modules<-list.files(my_target_dir)

dir.create('WGCNA_results/gprofiler_results')

save_fig_location<-'WGCNA_results/gprofiler_results/figures'
save_data_location<-'WGCNA_results/gprofiler_results/data_files'
dir.create(save_fig_location)
dir.create(save_data_location)

for (module in found_modules){
  
  module_name<-strsplit(module,'.csv')[[1]]
  save_name_fig<-paste0(save_fig_location,'/',module_name,'_overview.html')
  save_name_data<-paste0(save_data_location,'/',module_name,'_data.csv')
  
  gene_vect<-read.csv(paste0(my_target_dir,'/',module))
  gene_vect<-as.vector(gene_vect$geneSymbol)
  
  gostres <- gost(query = gene_vect,organism = "hsapiens")
  if (is.null(gostres)==F){
    #Save figure results
    p <- gostplot(gostres, capped = T, interactive = T)
    saveWidget(p,file=save_name_fig)
    
    #Produce csv results - conversion of parent column (collapse list)
    for_save<-gostres$result
    for (idx in 1:nrow(for_save)){
      for_save$parents[idx]<-paste(unlist(for_save$parents[idx]), collapse = '/')
    }
    for_save$parents<-unlist(for_save$parents)
    
    #Save data
    write.csv(x=for_save,file=save_name_data)
    
    #Delete extra folder created by htmlwidgets
    unlink(paste0(save_fig_location,'/',module_name,'_overview_files'),recursive = T)
  }
}

