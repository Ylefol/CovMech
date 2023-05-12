
setwd('~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/')
result_location<-'TS_results/gprofiler_results/data_files/'
list.files(result_location)

my_results<-list()
for(cluster in list.files(result_location)){
  cluster_id<-strsplit(cluster,'_data.csv')[[1]]
  temp_file<-read.csv(paste0(result_location,cluster),row.names = 1)
  temp_file<-temp_file[temp_file$source=='GO:BP',]
  temp_file<-temp_file[,c('source','term_name','term_id','p_value','term_size','query_size','intersection_size')]
  
  if(nrow(temp_file)>0){
    my_results[[cluster_id]]<-temp_file
  }
}

my_results<-my_results[c('C1','C5','C6','C7','C8','C9','C10')]

library(dplyr)
Cov_BP<-bind_rows(my_results, .id = "cluster_label")

write.csv(Cov_BP,'TS_cov_GO_BP.csv',row.names = F)
