
setwd('~/A_Projects/EpiGen/R_Work_Folder/Cov_Mech/')
library(ggraph)
library(ggplot2)
library(scales)
library(stringr)
library(DOSE)

load('data/GSEA_results.rds')

custom_dotplot_function <- function(gse_df, x = "geneRatio", color = "p.adjust",
                                    size=NULL, font.size=12, title = "", 
                                    orderBy="x",label_format = 50, decreasing=TRUE) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue","enrichmentScore"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    if (is.null(size))
      size  <- "Count"
  }
  

  #Replace activated and suppressed by better terminology
  gse_df$.sign<-gsub(x=gse_df$.sign,pattern = 'activated',replacement = 'positively enriched')
  gse_df$.sign<-gsub(x=gse_df$.sign,pattern = 'suppressed',replacement = 'negatively enriched')
  
  if (orderBy !=  'x' && !orderBy %in% colnames(gse_df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    gse_df <- dplyr::mutate(gse_df, x = eval(parse(text=x)))
  }
  
  idx <- order(gse_df[[orderBy]], decreasing = decreasing)

  new_col<-c()
  for(item in gse_df$Description){
    new_item<-paste(strwrap(item,60), collapse="\n")
    new_col<-c(new_col,new_item)
  }
  gse_df$Description<-new_col
  
  gse_df$Description <- factor(gse_df$Description,
                           levels=rev(unique(gse_df$Description[idx])))
  ggplot(gse_df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_fill_gradientn(colours=rev(c("red","white","blue")),
                         breaks=c(-1,0,1),
                         limits=c(-1, 1), oob=squish,
                         guide=guide_colorbar(reverse=F), name = color,aesthetics = "colour")+ 
    scale_y_discrete(labels = function(x) str_wrap(x, width = label_format))+
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
    scale_size(range=c(3, 8))
  
}

my_gg_save_function <- function(save_path_name,plot,height=NULL,width=NULL){
  tryCatch({
    if (is.null(height)==F){
      ggsave(save_path_name,dpi=300,height=height,width=width,unit='cm',
             plot,bg = "white")
    }else{
      ggsave(save_path_name,dpi=300,
             plot,bg = "white")
    }
    
  },
  error=function(cond){print(cond)}
  )
}


for(i in names(GSEA_res)){
  
  gse_df<-GSEA_res[[i]]
  
  gse_neg<-as.vector(gse_df$Description[gse_df$enrichmentScore<0][1:20])
  gse_pos<-as.vector(gse_df$Description[gse_df$enrichmentScore>0][1:20])
  gse_top_10_each<-c(gse_neg[1:10],gse_pos[1:10])  
  
  
  gse_pos_df<-gse_df[gse_df$Description %in% gse_pos,]
  my_gg_save_function(save_path_name = "dotplot_pos.png",width = 25, height = 20,
                      custom_dotplot_function(gse_pos_df, title = i,color="enrichmentScore") + facet_grid(.~.sign))
  
  gse_neg_df<-gse_df[gse_df$Description %in% gse_neg,]
  my_gg_save_function(save_path_name = "dotplot_neg.png",width = 25, height = 20,
                      custom_dotplot_function(gse_neg_df, title = i,color="enrichmentScore") + facet_grid(.~.sign))
  
  gse_both_df<-gse_df[gse_df$Description %in% gse_top_10_each,]
  my_gg_save_function(save_path_name = "dotplot_top10_each.png",width = 25, height = 20,
                      custom_dotplot_function(gse_both_df, title = i,color="enrichmentScore") + facet_grid(.~.sign))
  
  
  break #Break for test
}

#Data set up for merged dotplots
custom_GO_selection<-c('neutrophil mediated immunity','neutrophil degranulation','type I interferon signaling pathway',
                       'response to type I interferon','response to interferon-gamma','interferon-gamma-mediated signaling pathway',
                       'defense response to virus','cytokine-mediated signaling pathway','cytokine production',
                       'wound healing','coagulation','platelet degranulation','platelet activation','humoral immune response',
                       'RNA catabolic process','translation','ribosome biogenesis','antigen processing and presentation of peptide antigen via MHC class II',
                       'MHC class II protein complex assembly','T cell receptor signaling pathway')


#Generates different versions of the merged dotplots
target_exps<-c('all.cov.TP1_vs_healthy.controls','all.cov.TP2_vs_healthy.controls','all.cov.TP3_vs_healthy.controls')
target_exps<-c('severe.cov.TP1_vs_mild.cov.TP1','severe.cov.TP2_vs_mild.cov.TP2','severe.cov.TP3_vs_mild.cov.TP3')


names(target_exps)<-c('D1','D3','D8')
built_df<-data.frame(NULL)

for(i in 1:length(target_exps)){
  exp<-unname(target_exps)[i]
  timepoint<-names(target_exps)[i]
  
  temp_df<-GSEA_res[[exp]]
  temp_df<-temp_df[temp_df$Description %in% custom_GO_selection,c('ID','Description','setSize','enrichmentScore',"p.adjust")]
  temp_df[['-log10(padj)']]<--log10(temp_df$p.adjust)
  row.names(temp_df)=c(1:nrow(temp_df))
  temp_df$timepoint<-rep(x = timepoint,nrow(temp_df))
  
  if(nrow(built_df)==0){
    built_df<-temp_df
  }else{
    built_df<-rbind(built_df,temp_df)
  }
}
built_df$Description<-factor(built_df$Description,levels=rev(custom_GO_selection))
  
#Set-up necessary functions
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}


label_format<-20
label_func <- default_labeller(label_format)
if(is.function(label_format)) {
  label_func <- label_format
}


#With the full data, gotta make the plot now.
plt<-ggplot(built_df, aes(x=timepoint, y=Description, size=`-log10(padj)`, color=enrichmentScore)) +
  geom_point() +
  scale_fill_gradientn(colours=rev(c("red","white","blue")),
                       # values=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                       breaks=c(-1,0,1),
                       # breaks=c(1,0.80,0.60,0.40,0.20, 0,-0.20,-0.40,-0.60,-0.80,-1),
                       limits=c(-1, 1), oob=squish,
                       guide=guide_colorbar(reverse=F), name = 'enrichmentScore',aesthetics = "colour")+ 
  # scale_color_continuous(low="blue", high="gray", name = color,
  #                        limits=c(-1, 0),
  #                        guide='legend') +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 50)) +
  ylab(NULL) + ggtitle('All vs healthy') + theme_dose(15) +
  scale_size(range=c(3, 8))

ggsave(file='merged_tp_dotplot.png',plt,width = 10, height = 8,dpi=300)






