ml R/4.0.3
R
​
suppressMessages(library(data.table)) 
suppressMessages(library(tidyverse)) 
​
#Define trait info
dx_directory <- 'Cancer'
subfile_name <- 'Breast_2020'
#Define the filename the raw stats have
rawfile_here <- 'icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt'
trait_name <- gsub('.txt','',rawfile_here)
​
#Define the path where the raw sumstats live - change your paths to those on your personal computer
raw_path <- sprintf('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/resources/gwas_spredixcan_sumstats/%s/%s/gwas_raw/',dx_directory,subfile_name)
#Define where we want to save harmonized sumstats
outpath <- sprintf('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/resources/gwas_spredixcan_sumstats/%s/%s/gwas_for_spredixcan/',dx_directory,subfile_name)
​
#Read in stats
stats_raw <- as.data.frame(fread(paste0(raw_path,rawfile_here), quote = "",stringsAsFactors=F))
#Check header names 
head(stats_raw)
​
#Change chromosome from a numeric vector to a factor(category)
stats_raw$CHR <- as.factor(stats_raw$'chr.iCOGs')
stats_raw$POS <- as.numeric(stats_raw$'Position.iCOGs')
​
########################################
#Add GTEx rsIDs
hg_map <- data.frame(fread('/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/resources/GTEx_map_37.txt.gz',header=T,  sep="\t", quote = "",stringsAsFactors=F))
hg_map$CHR <- as.factor(hg_map$CHR)
hg_map$pos_37 <- as.numeric(hg_map$pos_37)
########################################
# #Filter and Join by chromosome and 37 position 
sumstats_merge_GTEx <- inner_join(stats_raw,hg_map,by=c('CHR','POS'='pos_37'))
head(sumstats_merge_GTEx)
​
#Remove variants without an rsID
sumstats_merge_GTEx_rsid <- sumstats_merge_GTEx %>% filter(!is.na(rs_id_dbSNP151_GRCh38p7))
########################################
#Document names of required s-predixcan columns in this data
#Rename pvalue without dots or slashes 
#Meta-analysis
​
sumstats_merge_GTEx_rsid$OR <- exp(sumstats_merge_GTEx_rsid$'Beta.meta')
​
sumstats_merge_GTEx_rsid <- sumstats_merge_GTEx_rsid %>% dplyr::rename(p='p.meta',se='sdE.meta')
​
#OR is included but no se
effect_size='OR'
effect_error='se'
​
########################################
#Determine effect allele and rename
​
sumstats_save_filt <- sumstats_merge_GTEx_rsid %>% dplyr::rename(effect_allele='Effect.Meta', noneffect_allele='Baseline.Meta')
​
# #Make sure reference/alternate alleles are the same between files
# # If GTEx snp is A/C we want the GWAS to be A/C not A/G
sumstats_save_filt <- sumstats_save_filt %>% filter((effect_allele == ref) | (effect_allele == alt)) %>% filter((noneffect_allele == ref) | (noneffect_allele == alt))
​
#Save new columns
sumstats_save_final <- sumstats_save_filt %>% select(all_of(c('rs_id_dbSNP151_GRCh38p7','effect_allele','noneffect_allele','p',effect_size,effect_error)))
​
#Combine into file name to later pull out in linux when we run s-predixcan
filename_new_here <- paste0(paste(subfile_name,trait_name,paste(effect_size,effect_error,sep='_'),sep='-'),'.txt.gz')
print(filename_new_here)
​
#Save
write.table(sumstats_save_final, file=gzfile(paste0(outpath,filename_new_here)), sep = "\t", col.names=T, row.names = F,quote=F)
