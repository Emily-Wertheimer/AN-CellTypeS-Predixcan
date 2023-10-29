## clear environment
rm(list=ls())
ls()

## install necessary packages if not already installed
if(!require(data.table)) {install.packages("data.table")}; library(data.table)
if(!require(tidyverse)) {install.packages("tidyverse")}; library(tidyverse)
if(!require(R.utils)) {install.packages("R.utils")}; library(R.utils)

## read in spredix results and tissues
spredix_software_path <-'~/AN-CellTypeS-Predixcan/sPrediXcan_AN/software'
path_predix_results_adiposeSub <-sprintf('~/AN-CellTypeS-Predixcan/sPrediXcan_AN/results/AN_GWAS_2019__Adipose_Subcutaneous.GTExv8')
predix_results_adiposeSub <- read.csv(path_predix_results_adiposeSub)
tissues <- read_lines('~/project/GTEx_v8_elastic_net_copies/en_Adipose_Subcutaneous.txt')
tissues_nospaces <- gsub(' ', '',tissues)
path_feature <- '/gpfs/gibbs/pi/huckins/projects/GTEx/BMI/featureData_bytissue/Adipose-Subcutaneous_featureData.txt.gz'
all_features <- NULL
tissue_here <- 'en_Adipose_Subcutaneous.txt'
path_feature_here <- sprintf('/gpfs/gibbs/pi/huckins/projects/GTEx/BMI/featureData_bytissue/Adipose-Subcutaneous_featureData.txt.gz')
featuredt <- as.data.frame(fread(path_feature_here, header=T, sep="\t", quote = "",stringsAsFactors=F))
featuredt$tissue <- tissue_here
all_features <- rbind(all_features,featuredt)

## left join predixcan w/ features to get gene start location
toplot_annot <- left_join(predix_results_adiposeSub,unique(all_features[,c('Name','start','Chr')]), by=c('gene' = 'Name'))

## add location col 
toplot_annot$location <- as.numeric(paste(toplot_annot$Chr, toplot_annot$start, sep = '.'))

## find bonf sig
bonf <- toplot_annot %>% dplyr::count() %>% mutate(bonf=.05/n) 
bonf_sig <- -log10(bonf$bonf)

# Create the Manhattan plot
plottingdir <- '/vast/palmer/home.mccleary/ekw28/AN-CellTypeS-Predixcan/sPrediXcan_AN/figures'
png("manhattan_plot.png", width = 12, height = 6, res = 300, units = "in")

manhattan_plot <- ggplot(toplot_annot, aes(x = as.numeric(location), y = -log10(pvalue), color=as.factor(Chr))) + 
  geom_point(size=2) + 
  theme_classic(base_size=24) +
  geom_hline(yintercept = bonf_sig, linetype='dashed', color='red') + 
  scale_x_continuous(limits=c(1,23),breaks=seq(1.5,22.5), labels=seq(1,22)) + 
  labs(x='Chromosome', y='-log10(p)')
dev.off()



# Display the plot
print(manhattan_plot)

