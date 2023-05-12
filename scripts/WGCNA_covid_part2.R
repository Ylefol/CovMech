options(stringsAsFactors = FALSE)
setwd("~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/")
source("R/my_WGCNA_functions.R")


target_power=18

#Still gotta load the R data
load("WGCNA_results/WGCNA_part1.RData")

# 
#The below function has many parameters and adapting them to our data
#May be necessary. Especially the maxBlockSize as this can 'break'
#the downstream code if not calibrated properly.
net = blockwiseModules(datExpr, power = target_power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       maxBlockSize = 35000,
                       verbose = 3,nThreads=4)

# saveTOMFileBase = paste0("WGCNA_results/",target_exp,"/block_wise_modules"),



save.image("WGCNA_results/WGCNA_block_calc_temp.RData")

load("WGCNA_results/WGCNA_block_calc_temp.RData")

# We now display the hierachical dendrograme
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath

png("WGCNA_results/dendrograme_plus_colors.png",750,500)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# We note that the function called 'recutBlockwiseTrees' can allow us to modify
# the dendrograme without having to re-compute everything.


#Renaming for downstream code
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;

# Define numbers of genes and samples
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#Using trycatch to circumvent an error with ggsave which does not impact the 
#figure but will halt the script. Using ggsave instead of png to easily
#obtain a higher resolution for this plot as it is quite packed with information.
png("WGCNA_results/labelled_heat_map.png",25,15,units='cm',res=600)
plot_labeled_heatmap(moduleTraitCor,moduleTraitPvalue,datTraits,MEs)
dev.off()
# plot_labeled_heatmap(net,datExpr,datTraits)

gene_names=colnames(datExpr)

dir.create("WGCNA_results/coexpression_genes")


# for (i in colnames(datTraits)){
#   csv_res_single_clinical(clin_dat = i,datTraits = datTraits,datExpr = datExpr,
#                           MEs=MEs,nSamples=nSamples,moduleColors=moduleColors,
#                           gene_names=gene_names)
#   
# }

dir.create('WGCNA_results/WGCNA_module_gene_lists/')
module_gene_info = data.frame(geneSymbol = gene_names, moduleColor = moduleColors)
modules<-unique(module_gene_info$moduleColor)
for (module in modules){
  extracted_genes<-module_gene_info[module_gene_info$moduleColor==module,]
  write.csv(extracted_genes,paste0('WGCNA_results/WGCNA_module_gene_lists/',module,'.csv'))
}


save.image("WGCNA_results/WGCNA_part2.RData")#Saves all objects in memory in a binary file
  



