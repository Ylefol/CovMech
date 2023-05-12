# CovMech
## Data, functions, and scripts for the CoVMech (Covid mechanism) project

The repository contains a variety of scripts which are executed in different ways.\

The first is a time series analysis performed with the 'TS_CoV_script.Rmd' which utilizes Rmarkdown to run an analysis. The Rmarkdown file itself explains the way to utilize it. This TimeSeries analysis is an in-development version of what has now become [TiSA](https://github.com/Ylefol/TimeSeriesAnalysis) which has also been [published](https://academic.oup.com/nargab/article/5/1/lqad020/7069286). The code is presented in this repository as the published version has since been changed since the analysis of this data. However if one were to use the time series pipeline, utilizing the published version is recommended as the published version is better adapted to various datasets and contains less bugs. It is also being maintained while the code here will not be maintained. The code in this repository exists for the exclusive purpose of sharing the code for the manuscript relating to the covid analysis.\

#### scripts folder\
- **ClusterProfiler_calls.R** - Performs a GSEA analysis\
- **Limma_processing.R** - Performs processing of data with Limma, the results of this processing are already located in the 'data' and therefore it does not need to be re-run. The reason for this is that the inclusion of each microarray file would results in a ~3GB repository. The files can be found with the following GSE identifier: GSE213313\
- **Limma_DE_cov.R** - Performs differential gene expression with limma and produces a set of PCA plots per comparison\
- **WGCNA_covid_part1.R** - Performs the first step of a WGCNA analysis\
- **WGCNA_covid_part2.R** - performs the second step of a WGCNA analysis\

#### script sub-folder\
- **clin_dta_traj_scripts.R** - This script creates trajectory plots for the given clinical parameters\
- **cov_extract_GO_BP.R** - Creates a csv file containing the GOs found within clusters identified by the time series analysis\
- **covid_article_plots.R** - Generates several plots, primarily it creates dotplots of custom GOs shown at each requested comparison. It also creates a heatmap for the time series analysis as well as a excel file with the GOs found in the time series analysis.\
- **covid_clust_traj_healthy.R** - Creates trajectory plots of either clusters of time series analysis or WGCNA modules. It also adds green lines representing healthy patients.\
- **covid_health_heat.R** - Creates a custom heatmap comparing covid patients with healthy controls.\
- **CoV_specific_figs.R** - Creates a correlation heatmap for WGCNA data and dotplots for GSEA data. The dotplots show the top 10 positively enriched and top 10 negatively enriched GOs of each comparison provided.\
- **create_dotplots_from_data.R** - This folder creates the dotplots from pre-saved data. This script will reproduce the dotplots as seen in the manuscript. The reason that the pre-saved data needs to be used is that updates to the GO databases has made it so that re-runs marginally changes the results. Overall interpretation would not change with the updated database, however small visual difference would be seen.\
- **WGCNA_gprofiler.R** - Performs gprofiler analysis for the identified modules within WGCNA\


As a generality, the scripts in the subfolder require the results of DEA, GSEA, WGCNA, and time series, therefore these should be run prior to the sub-folder scripts.
