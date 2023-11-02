## load necessary libraries
if(!require(data.table)) {install.packages("data.table")}; library(data.table)
if(!require(tidyverse)) {install.packages("tidyverse")}; library(tidyverse)
if(!require(R.utils)) {install.packages("R.utils")}; library(R.utils)

## load data
stats_raw <- as.data.frame(fread('/home/ekw28/AN-CellTypeS-Predixcan/sPrediXcan_AN/data/freeze3sumStats/daner_AN.meta.gz', quote="", stringsAsFactors=F)) 

hg_map <- data.frame(fread('/gpfs/gibbs/pi/huckins/resources/GTEx_map_37.txt.gz'
,header=T,  sep="\t", quote = "",stringsAsFactors=F))

## change chr & pos to factors
stats_raw$CHR <- as.factor(stats_raw$CHR) 
stats_raw$BP <- as.numeric(stats_raw$BP)

hg_map$CHR <- as.factor(hg_map$CHR) 

## id & rename ref & alt 
colnames(stats_raw)[4]<-'ref'
colnames(stats_raw)[5]<-'alt'

## join sumstats w/ human genome map by chromosome & position 
sumstats_merge_GTEx <- inner_join(stats_raw,hg_map,by=c('CHR','BP'='pos_37'))

colnames(sumstats_merge_GTEx_rsid)
# [1] "CHR"                     "SNP"                    
 #[3] "BP"                      "ref.x"                  
# [5] "alt.x"                   "FRQ_A_24145"            
 #[7] "FRQ_U_1244243"           "INFO"                   
 #[9] "OR"                      "SE"                     
#[11] "P"                       "ngt"                    
#[13] "Direction"               "HetISqt"                
#[15] "HetDf"                   "HetPVa"                 
#[17] "Nca"                     "Nco"                    
#[19] "Neff_half"               "ref.y"                  
#[21] "alt.y"                   "variant_id_hg38"        
#[23] "rs_id_dbSNP151_GRCh38p7"


colnames(sumstats_merge_GTEx_rsid)[4] <- 'noneffect_allele'
colnames(sumstats_merge_GTEx_rsid)[5] <- 'effect_allele'
colnames(sumstats_merge_GTEx_rsid)[11] <- 'p'
colnames(sumstats_merge_GTEx_rsid)[22] <- 'variant.id.hg38'
colnames(sumstats_merge_GTEx_rsid)[23] <- 'rs.id.dbSNP151.GRCh38p7'
colnames(sumstats_merge_GTEx_rsid)[10] <- 'effect_error'

dim(sumstats_merge_GTEx_rsid) #dim: 6804937      23


######################################## 

sumstats_save_filt <- sumstats_merge_GTEx_rsid # %>% dplyr::rename(effect_allele='alt', noneffect_allele='ref')

#Save new columns
sumstats_save_final <- sumstats_save_filt %>% select(all_of(c('rs.id.dbSNP151.GRCh38p7','effect_allele','noneffect_allele','p','OR')))

#Combine into file name to later pull out in linux when we run s-predixcan
#filename_new_here <- paste0(paste(subfile_name,trait_name,paste(effect_size,effect_error,sep='_'),sep='-'),'.txt.gz')
#print(filename_new_here)
#couldnt get previous 2 lines to work so did this instead: 
freeze3.gwas.sumstats.clean <- sumstats_save_final
​head(anGWAS2019.sumstats.clean)

#Save
#write.table(freeze3.gwas.sumstats.clean, file = "/vast/palmer/home.mccleary/ekw28/freeze3.gwas.sumstats.clean")
write.table(freeze3.gwas.sumstats.clean, file=/home/ekw28/AN-CellTypeS-Predixcan/freeze3.gwas.sumstats.clean2, sep = "\t", col.names=T, row.names = F,quote=F)
