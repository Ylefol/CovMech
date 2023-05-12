setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")

clin_dta<-read.csv('data/clinical_data_for_traj.csv',row.name=1)


#Convert to numeric and log transform
clin_dta$RNA <- gsub(pattern = ',',replacement = '.',x = clin_dta$RNA)
clin_dta$RNA <- as.numeric(clin_dta$RNA)
clin_dta$RNA<-log10(clin_dta$RNA)
clin_dta$RNA <- gsub(pattern = '-Inf',replacement = 0,x = clin_dta$RNA)

clin_dta$IgG <- as.numeric(clin_dta$IgG)
clin_dta$IgG<-log10(clin_dta$IgG)
clin_dta$IgG <- gsub(pattern = '-Inf',replacement = 0,x = clin_dta$IgG)

clin_dta$totIg<- as.numeric(clin_dta$totIg)
clin_dta$totIg<-log10(clin_dta$totIg)
clin_dta$totIg <- gsub(pattern = '-Inf',replacement = 0,x = clin_dta$totIg)

clin_dta$Neutrophils<- as.numeric(clin_dta$Neutrophils)
clin_dta$Neutrophils<-log10(clin_dta$Neutrophils)
clin_dta$Neutrophils <- gsub(pattern = '-Inf',replacement = 0,x = clin_dta$Neutrophils)

library(reshape2)
clin_dta$ID<-unlist(strsplit(clin_dta$ID,'_'))[c(T,F)]#Split and keep uneven (ID, not _TP)

melted_vals <- reshape2::melt(clin_dta,id=c('ID','Timepoint','Severe'))
melted_vals <- melted_vals[is.na(melted_vals$value)==F,]#Remove NAs to prevent plotting
melted_vals$value<-as.numeric(melted_vals$value)
# test$value<-round(as.numeric(test$value), digits = 2)
# test<-test[test$variable != 'ID',]


graphic_vector<-c("#e31a1c","#1f78b4")
names(graphic_vector)<-c('Critical','Non-critical')

melted_vals$Severe[melted_vals$Severe=='Y']<-"Critical"
melted_vals$Severe[melted_vals$Severe=='N']<-"Non-critical"

names(melted_vals)[names(melted_vals) == 'Severe'] <- 'Group'


#Convert timepoints to numeric
melted_vals$Timepoint[melted_vals$Timepoint=='T1']<-0
melted_vals$Timepoint[melted_vals$Timepoint=='T2']<-48
melted_vals$Timepoint[melted_vals$Timepoint=='T3']<-168


melted_vals$Timepoint<-factor(as.numeric(melted_vals$Timepoint),levels=c(0,48,168))
#Create custom labels
melted_vals$labels<-paste0(melted_vals$variable,' – ',melted_vals$Group)

#Add number of patients at each timepoint to the labels
for (label in unique(melted_vals$labels)){
  num_samples<-table(melted_vals$Timepoint[melted_vals$labels==label])
  new_label<-paste0(label,': Day1 – ',unname(num_samples[1]),' | Day3 – ',unname(num_samples[2]),' | Day8 – ',unname(num_samples[3]))
  melted_vals$labels[melted_vals$labels==label]<-new_label
}


library(ggplot2)

for(dta_type in unique(melted_vals$variable)){
  
  sub_melt<-melted_vals[melted_vals$variable == dta_type,]
  
  #calculate mean
  mean_df<-data.frame(NULL)
  for(cat in unique(sub_melt$labels)){
    temp_df<-sub_melt[sub_melt$labels==cat,]
    for(tp in unique(temp_df$Timepoint)){
      mean_val<-mean(temp_df$value[temp_df$Timepoint==tp])
      df<-data.frame(Timepoint=tp,Group=unique(temp_df$Group),value=mean_val,labels=cat)
      if(nrow(mean_df)==0){
        mean_df<-df
      }else{
        mean_df<-rbind(mean_df,df)
      }
    }
  }
  
  #Conditional for renaming Y axis
  if(dta_type=='IgG'){
    y_axis_name<-'Anti-SARS-CoV-2 Spike IgG \n (log10 AU/mL)'
  }else if(dta_type=='totIg'){
    y_axis_name<-'Total anti-SARS-CoV-2 NP Ig \n (log10 AU/mL)'
  }else{
    y_axis_name<-dta_type
  }
  
  
  
  plt <- ggplot(sub_melt, aes(y = value , x =Timepoint, color = Group))+
    scale_color_manual(values=graphic_vector) +
    geom_line(aes(group = ID), alpha = 0.4) +
    geom_point() +
    ylab('log10(values)')+
    geom_line(
      data = mean_df, lwd = 1.5, color = "grey50",
      aes(group = labels)
    ) +
    scale_x_discrete(name ="timepoints",
                     breaks=c('0','48','168'),
                     labels=c("D1", "D3", "D8"))+
    ylab(y_axis_name)+
    facet_wrap(~labels,ncol = 2)
  
  svg(paste0(dta_type,'.svg'),width=8,height=3)
  print(plt)
  dev.off()
  
  
}


