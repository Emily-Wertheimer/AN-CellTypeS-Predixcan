# ## clear environment
rm(list=ls())
ls()
# 
## install necessary packages if not already installed
if(!require(data.table)) {install.packages("data.table")}; library(data.table)
if(!require(tidyverse)) {install.packages("tidyverse")}; library(tidyverse)
if(!require(R.utils)) {install.packages("R.utils")}; library(R.utils)
if(!require(ggplot2)) {install.packages("ggplot2")}; library(ggplot2)

## define tissue
tissue_predix <- 'dACC'
cell_type_predix <- '_Excitatory_Neurons'
suffix <- '.carina.cell.type' 

## read in spredix results and tissues
spredix_software_path <-'/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py' # ~/AN-CellTypeS-Prredixcan/sPrediXcan_AN/software'
path_spredix_results <- sprintf(paste0('/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/', tissue_predix, cell_type_predix, suffix))
predix_results <- read.csv(path_spredix_results)
tissues <- read_lines('/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/Excitatory_Neurons/dACC_Excitatory_Neurons_LD_covariance_matrix.txt.gz')
tissues_nospaces <- gsub(' ', '',tissues)
 
all_features <- NULL
tissue_here <- 'dACC_Excitatory_Neurons'
path_feature_here <- sprintf('/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/carina_cellType_annotations.txt')# featuredt <- as.data.frame(fread(path_feature_here))#, header=T, sep="\t", quote = "",stringsAsFactors=F))
featuredt <- as.data.frame(fread(path_feature_here))
featuredt$tissue <- tissue_here
all_features <- rbind(all_features,featuredt)

## left join predixcan w/ features to get gene start location
toplot_annot <- left_join(predix_results,unique(all_features[,c('Gene_name','Gene_start','Chromosome')]), by=c('gene_name'='Gene_name'))

## add location col 


# export
toplot_annot$location <- gsub("[^0-9]", "", toplot_annot$Chromosome) #removes "chr" (e.g. "chr3" --> "3")
toplot_annot$locationBP <- as.numeric(paste(toplot_annot$location, toplot_annot$Gene_start, sep='.'))

## find bonf sig
bonf <- toplot_annot %>% dplyr::count() %>% mutate(bonf=.05/n) 
bonf_sig <- -log10(bonf$bonf)

# Create the Manhattan plot
plottingdir <- '/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/figures'
png("iGlut_dACC_manhattan_plot.png", width = 12, height = 6, res = 300, units = "in")

manhattan_plot <- ggplot(toplot_annot, aes(x = as.numeric(locationBP), y = -log10(pvalue), color=as.factor(Chromosome))) + 
  geom_point(size=2) + 
  theme_classic(base_size=24) +
  geom_hline(yintercept = bonf_sig, linetype='dashed', color='red') + 
  scale_x_continuous(limits=c(1,23),breaks=seq(1.5,22.5), labels=seq(1,22)) + 
  labs(x='Chromosome', y='-log10(p)')
dev.off()



# Display the plot
print(manhattan_plot)
