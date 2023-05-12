
setwd(dir = "~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")
source("R/my_ClusterProfiler_functions.R")

lst<-list.files('DE_results/')
lst<-lst[! lst %in% c('DE_summary.csv','Robjects_checkpoints')]
# file_to_use <- 'DE_results_with_gtf.csv'
file_to_use <- 'DE_raw_data.csv'

library(org.Hs.eg.db)
selected_org_GO="org.Hs.eg.db"
selected_org_gsea=org.Hs.eg.db


my_filter='padj'
clustpro_path<-'cluster_profiler_results'
de_path<-'DE_results'
#Create the results folder
dir.create(clustpro_path)


#Set up the folders and call the wrapper analysis
library(GOSemSim)
my_go_dta <- godata(selected_org_GO, ont="BP")
for (item in lst){
  if (item %in% list.files(clustpro_path)==TRUE){
    next
  }
  DE_file <- read.csv(paste0(de_path,"/",item,"/",file_to_use),check.names = FALSE)
  dir.create(paste0(clustpro_path,"/",item))
  dir.create(paste0(clustpro_path,"/",item,"/GSEA"))
  save_path=paste0(clustpro_path,"/",item)

  wrapper_for_GSEA_and_overrepresentation(DE_file,filter_choice=my_filter,save_path,selected_org_gsea,selected_org_GO,go_dta=my_go_dta,exp_name=item)
}

