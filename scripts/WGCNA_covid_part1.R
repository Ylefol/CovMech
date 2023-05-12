
#Usefull ressource if I encounter problems
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html#:~:text=How%20many%20samples%20do%20I,network%20to%20be%20biologically%20meaningful.

library('WGCNA')
library(org.Hs.eg.db)

options(stringsAsFactors = FALSE)
setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech//")
source("R/my_WGCNA_functions.R")

dir.create("WGCNA_results")

#Load data
datExpr0<-read.csv('data/WGCNA_dta/WGCNA_datExpr0.csv',row.names = 1)
datTraits<-read.csv('data/WGCNA_dta/WGCNA_datTraits.csv',row.names = 1)

#Do this so that the row order matches between them
datExpr0 <- datExpr0[sort(rownames(datExpr0)),]
datTraits <- datTraits[sort(rownames(datTraits)),]

# We then do some QC by checking for outliers and removing samples with too many missing values
datExpr0 <- pre_WGCNA_QC(datExpr0)

# We now cluster the samples
sampleTree = hclust(dist(datExpr0), method = "average")


#The plot shows potential outliers, we could remove them by hand or automatically
#The preferred solution is to do so automatically and check visually if any remain after
#having implemented the proper solution
png("WGCNA_results/Sample_tree.png",width = 800 ,height=720)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

# Determine cluster under the line, re-cluster if necessary
returned_list <- cluster_under_line(sample_tree,datExpr0,datTraits)

#Assign the returned values to variables
sampleTree2 <- returned_list[[1]]
traitColors <- returned_list[[2]]
datExpr <- returned_list[[3]]

#Save the clinical data with associated dendrogram
svg("WGCNA_results/clinical_with_dendrograme.svg",15,10)
# Plot the sample dendrogram and the colors underneath.
print(plotDendroAndColors(sampleTree2, traitColors,
                          groupLabels = names(datTraits),
                          main = "Sample dendrogram and trait heatmap"))
dev.off()

# source("code/functions/my_WGCNA_functions.R")
svg(paste0("WGCNA_results/clinical_values_heatmap.svg"),15,10)
print(plot_clinical_heatmap(datTraits,traitColors,'Critical'))
dev.off()


#Save the power plots for subsequent user analysis
png(paste0("WGCNA_results/power_plot.png"),width=750,height=750)
print(plot_power(datExpr))

dev.off()

library(corrplot)
# setwd("~/A_Projects/EpiGen/R_Work_Folder/RNAseq_pipeline/")
# load("WGCNA_results/WGCNA_part1.RData")
# source('code/functions/my_WGCNA_functions.R')


p.mats<-calculate_corr_and_pval_matrix(datTraits)
my_format_p<-format_p_val_matrix(p.mats[[1]])
corr_mat<-p.mats[[2]]
plot_clinical_corr_plot(datTraits,corr_mat,my_format_p)

my_format_p[my_format_p != ""]<-'*'
plot_clinical_corr_plot_V2(datTraits,corr_mat,my_format_p)

#Save the object from this script so it may be taken by the second step of the WGCNA analysis
save.image(paste0("WGCNA_results/WGCNA_part1.RData"))#Saves all objects in memory in a binary file

       