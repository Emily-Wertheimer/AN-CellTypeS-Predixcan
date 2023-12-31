## get started
#Download AN 2019 GWAS sumstats from pgc website (don't save!)
#Upload to McCleary (gpfs/gibbs/project/levy_ifat/ekw28/) using YCRC OOD GUI
#request resources
#load R

##Install packagaes 
install.packages('data.table')
library('data.table')
install.packages('dplyr')
library(dplyr)
install.packages('tidyverse')
library(tidyverse)

########################################
#Define the path where the raw sumstats live - change your paths to those on your personal computer
raw_path <- sprintf('/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/data/anGWAS2019.download.orig.txt')

#Define where we want to save harmonized sumstats
outpath <- sprintf('/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan.AN/data/anGWAS2019.sumstats.clean.2')

#Read in stats
stats_raw <- as.data.frame(fread('/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/data/anGWAS2019.download.orig.txt', quote="", stringsAsFactors=F)) #dim: 8219102      14

#Change chromosome from a numeric vector to a factor(category)
stats_raw$CHROM <- as.factor(stats_raw$CHROM) 
stats_raw$POS <- as.numeric(stats_raw$POS) 

########################################
#Add GTEx rsIDs
hg_map <- data.frame(fread('/gpfs/gibbs/pi/huckins/resources/GTEx_map_37.txt.gz',header=T,  sep="\t", quote = "",stringsAsFactors=F)) #dim 465697 6
head(hg_map)
hg_map$CHR <- as.factor(hg_map$CHR) 
hg_map$pos_37 <- as.numeric(hg_map$pos_37) 

########################################
#Filter and Join by chromosome and 37 position 
colnames(stats_raw)[1] <- 'CHR' # make sure colnames match bw two matrices, #dim: 1048505      14
stats_raw$CHR <- as.factor(stats_raw$CHR) # can check using class() #dim: 104805 14
sumstats_merge_GTEx <- inner_join(stats_raw,hg_map,by=c('CHR','POS'='pos_37'))
head(sumstats_merge_GTEx)
# warning message: In inner_join(stats_raw, hg_map, by = c("CHR", POS = "pos_37")) :
 # Detected an unexpected many-to-many relationship between `x` and `y`.
 #ℹ Row 383 of `x` matches multiple rows in `y`.
 #ℹ Row 42548644 of `y` matches multiple rows in `x`.
 #ℹ If a many-to-many relationship is expected, set `relationship = "many-to-many"` to silence this warning.
head(sumstats_merge_GTEx)​
dim(sumstats_merge_GTEx) #dim: 7422616      18


#Remove variants without an rsID
sumstats_merge_GTEx_rsid <- sumstats_merge_GTEx %>% filter(!is.na(rs_id_dbSNP151_GRCh38p7)) # doesn't look like this filtered anything

########################################
#Document names of required s-predixcan columns in this data
  # rs_id_dbSNP151_GRCh38p7
  # effect_allele
  # noneffect_allele
  # p 
  # OR (instead of beta)

##Determine effect allele (from GWAS readme, confirm using paper fig) and rename
  # 2019 AN GWAS: effect is given in beta. exp(beta) = odds ratio 
  # ref = A2
  # alt = A1
  # odds are in setting of alt allele

#Rename pvalue without dots or slashes 
colnames(sumstats_merge_GTEx_rsid)
# [1] "CHR"                     "POS"                    
# [3] "ID"                      "REF"                    
# [5] "ALT"                     "BETA"                   
# [7] "SE"                      "PVAL"                   
# [9] "NGT"                     "IMPINFO"                
# [11] "NEFFDIV2"                "NCAS"                   
# [13] "NCON"                    "DIRE"                   
# [15] "ref"                     "alt"                    
# [17] "variant_id_hg38"         "rs_id_dbSNP151_GRCh38p7"
colnames(sumstats_merge_GTEx_rsid)[4] <- 'noneffect_allele'
colnames(sumstats_merge_GTEx_rsid)[5] <- 'effect_allele'
colnames(sumstats_merge_GTEx_rsid)[8] <- 'p'
colnames(sumstats_merge_GTEx_rsid)[17] <- 'variant.id.hg38'
colnames(sumstats_merge_GTEx_rsid)[18] <- 'rs.id.dbSNP151.GRCh38p7'
colnames(sumstats_merge_GTEx_rsid)[7] <- 'effect_error'
colnames(sumstats_merge_GTEx_rsid)[19] <- 'OR'

dim(sumstats_merge_GTEx_rsid) #dim: 945197     18

# change BETA to OR and rename
head(sumstats_merge_GTEx_rsid$BETA)
# [1]  0.009197572  0.000199980  0.010900374  0.027595712 -0.024302939
# [6]  0.018399683
sumstats_merge_GTEx_rsid$OR <- exp(sumstats_merge_GTEx_rsid$'BETA')
head(sumstats_merge_GTEx_rsid$OR)
# [1] 1.00924 1.00020 1.01096 1.02798 0.97599 1.01857

# check odds calculation
#pick SNP from table 1 in 2019 paper: SNP = rs9821797, OR = 1.17, s.e. = 0.02
filter(sumstats_merge_GTEx_rsid, ID=='rs62513865') # filter for chosen snp
exp(#beta of chosen snp)
#dim: 945197 19
#colnames
# [1] "CHR"                     "POS"                    
# [3] "ID"                      "noneffect_allele"       
# [5] "effect_allele"           "BETA"                   
# [7] "SE"                      "p"                      
# [9] "NGT"                     "IMPINFO"                
# [11] "NEFFDIV2"                "NCAS"                   
# [13] "NCON"                    "DIRE"                   
# [15] "ref"                     "alt"                    
# [17] "variant.id.hg38"         "rs.id.dbSNP151.GRCh38p7"
# [19] "OR"   

######################################## 

sumstats_save_filt <- sumstats_merge_GTEx_rsid # %>% dplyr::rename(effect_allele='alt', noneffect_allele='ref')

# #Make sure reference/alternate alleles are the same between files
# # If GTEx snp is A/C we want the GWAS to be A/C not A/G
sumstats_save_filt <- sumstats_save_filt %>% filter((effect_allele == ref) | (effect_allele == alt)) %>% filter((noneffect_allele == ref) | (noneffect_allele == alt))
dim(sumstats_save_filt) # 6911447      19

#Save new columns
sumstats_save_final <- sumstats_save_filt %>% select(all_of(c('rs.id.dbSNP151.GRCh38p7','effect_allele','noneffect_allele','p','OR')))
​ # effect size = OR

#Combine into file name to later pull out in linux when we run s-predixcan
#filename_new_here <- paste0(paste(subfile_name,trait_name,paste(effect_size,effect_error,sep='_'),sep='-'),'.txt.gz')
#print(filename_new_here)
#couldnt get previous 2 lines to work so did this instead: 
anGWAS2019.sumstats.clean.2 <- sumstats_save_final
​head(anGWAS2019.sumstats.clean)

#Save
outpath <- sprintf('/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/data/anGWAS2019.sumstats.clean.2.txt')
write.table(anGWAS2019.sumstats.clean.2, file=outpath, sep = "\t", col.names=T, row.names = F,quote=F)
