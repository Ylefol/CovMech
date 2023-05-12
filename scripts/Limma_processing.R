library(limma)
library(HsAgilentDesign026652.db)
setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")

# read targets and data files
targets<-readTargets("data/micro_arr_metadata.txt", row.names="ID")


# NMLR <- sqrt(targets$Neutrophils*targets$Monocytes)/targets$Lymphocytes
# targets$NMLR <- NMLR
targets$Severe[targets$Group == "Ctrl"] <- "C"
targets$MechVent[targets$Group == "Ctrl"] <- "C"
targets$Viremi[targets$Group == "Ctrl"] <- "C"
targets$Sex[targets$Sex == 1] <- "M"
targets$Sex[targets$Sex == 0] <- "F"

setwd('data/Microarray_data')
arraydata <- read.maimages(targets, source="agilent", names=targets$ID, green.only=TRUE, other.columns=c("gIsWellAboveBG", "gIsFeatNonUnifOL", "gIsBGNonUnifOL", "gIsFeatPopnOL", "gIsBGPopnOL"))

# add annotations
arraydata$genes$EntrezID <- mapIds(HsAgilentDesign026652.db, arraydata$genes$ProbeName, keytype="PROBEID", column="ENTREZID")
# arraydata$genes$Symbol   <- mapIds(HsAgilentDesign026652.db, arraydata$genes$ProbeName, keytype="PROBEID", column="SYMBOL")

# use normexp background correction followed by quantile normalization:
arraydata.normalized <- backgroundCorrect(arraydata, method="normexp")
arraydata.normalized <- normalizeBetweenArrays(arraydata.normalized, method="quantile")

# save("arraydata.normalized", file="my_covid_data_normalized.rdata")

# filter out control probes, outliers, and probes that don’t appear to be expressed.
# keep probes that are above background or not outliers at least 1/4 of the arrays:
Control <- arraydata.normalized$genes$ControlType==1L
IsExpr <- rowSums(arraydata.normalized$other$gIsWellAboveBG > 0) >= ncol(arraydata.normalized)/4
isOutlier <- rowSums(arraydata.normalized$other$gIsFeatNonUnifOL == 1 |
                       arraydata.normalized$other$gIsBGNonUnifOL == 1 |
                       arraydata.normalized$other$gIsFeatPopnOL == 1 |
                       arraydata.normalized$other$gIsBGPopnOL == 1) >= ncol(arraydata.normalized)/4

# select the probes to keep in a new data object:
arraydata.filtered <- arraydata.normalized[!Control & IsExpr &  !isOutlier, ]

# Average duplicate probes
arraydata.filtered <- avereps(arraydata.filtered, ID=arraydata.filtered$genes$SystematicName)
#arraydata.filtered <- avereps(arraydata.filtered, ID=arraydata.filtered$genes$GeneName)

# remove annotation columns we no longer need:
arraydata.filtered$genes <- arraydata.filtered$genes[,c("ProbeName","GeneName","SystematicName","EntrezID","Description")]

#Make submission matrix
# gene_df<-data.frame(ID_REF=arraydata.filtered$genes$ProbeName)
# 
# val_df<-data.frame(arraydata.filtered$E)
# colnames(val_df)=colnames(arraydata.filtered$E)

# final_sub_df<-cbind(gene_df,val_df)
# library("xlsx")
# openxlsx::write.xlsx(final_sub_df, 'Metadata_Matrix.xlsx', sheetName = "Matrix",
#                      colNames = TRUE, rowNames = FALSE, append = FALSE)

library(dplyr)
#Rename
targets %>%
  dplyr::rename(
    SampleID = ID,
  ) -> targets

setwd('~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/')
rna_biop_dataa<-arraydata.filtered
patient_sample_dta<-targets
save(rna_biop_dataa,patient_sample_dta, file='data/processed_microarray_data.rdata', version = 2)






