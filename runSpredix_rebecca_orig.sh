#Run s-predixcan on ANM and BEB CCGWAS - exact pvalue

salloc -t 6:00:00 --mem=32G 
#r206u22n04
#Local R library is: /home/rs2826/project/R

ml R/4.2.0-foss-2020b

cat("SETUP: loading libraries\n")
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

#Read in results.all file
stat_raw <- as.data.frame(fread('/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/results/ANM_BEB_ccgwas.results.ALL.gz', header=T, sep="\t", quote = "",stringsAsFactors=F))
head(stat_raw)

######################################################
#Prep file for s-predixcan
stat_rename <- stat_raw %>% dplyr::rename(effect_allele=EA,noneffect_allele=NEA,p=Exact_pval, beta=Exact_beta) %>% mutate(CHR=as.factor(CHR)) 

######################################################
#Merge to GTEX v8 rsIDs
#Read in a mapping of hg37 positions to hg38
hg37_map <- data.frame(fread('/gpfs/gibbs/pi/huckins/resources/GTEx_map_37.txt.gz',header=T,  sep="\t", quote = "",stringsAsFactors=F))
#Change chromosome from a numeric vector to a factor(category)
hg37_map$CHR <- as.factor(hg37_map$CHR)

#Filter and Join by chromosome and 37 position 
sumstats_merge_GTEx <- inner_join(stat_rename,hg37_map,by=c('CHR','BP'='pos_37'))

#Make sure reference/alternate alleles are the same between files
# If GTEx snp is A/C we want the GWAS to be A/C not A/G
sumstats_merge_GTEx_refalt <- sumstats_merge_GTEx %>% filter((effect_allele == ref) | (effect_allele == alt)) %>% filter((noneffect_allele == ref) | (noneffect_allele == alt))

#Remove variants without an rsID
sumstats_merge_GTEx_rsid <- sumstats_merge_GTEx_refalt %>% filter(!is.na(rs_id_dbSNP151_GRCh38p7))

final <- sumstats_merge_GTEx_rsid %>% select(rs_id_dbSNP151_GRCh38p7, CHR, effect_allele, noneffect_allele, beta,p)

#Save stats
outpath_spredix <- '/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/results/ANM_BEB_ccgwas.results.ALL_prepped_spredix_exact.gz'
write.table(final, file=gzfile(outpath_spredix),quote=F, sep="\t", col.names=T, row.names=F)

#Quit

##############################################################################
##############################################################################

#Create conda environment

# module purge && module load miniconda
# conda env create -n py3_env_spredix --file /gpfs/gibbs/pi/huckins/software/MetaXcan/software/conda_env.yaml
# module purge && module load miniconda
# conda activate py3_env_spredix

#Run job 

sbatch /gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/jobs/GTExv8_loop_ANM_BEB_exact_rsmap.sh

####################################################################
####################################################################

#!/bin/bash
#SBATCH --job-name=GTExv8_loop_ANM_BEB_exact_rsmap
#SBATCH --output=/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/scratch/%x
#SBATCH --time=4:00:00
#SBATCH --mem 64G 

module purge && module load miniconda
conda activate py3_env_spredix

#Run s-predixcan on ccgwas output
SPREDIX_path=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py
#Define outpath and scratch
outpath=/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/spredixcan
scratch=/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/scratch

#Trait and GWAS file name and location
trait=ANM_BEB
p_type=exact
path_GWAS=/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/results/
trait_filename=${trait}_ccgwas.results.ALL_prepped_spredix_${p_type}.gz

#Loop and run
LOOPINGFILE=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/Tissues
cd $scratch
while read tissue_predix; do
#Define tissue
#tissue_predix=Thyroid
echo $tissue_predix

#Define location of models
PREDICT_DB_path=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_${tissue_predix}.db
COVARIANCE_PATH=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_${tissue_predix}.txt.gz
jobname=${trait}_${tissue_predix}_spredixcan

python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename} --snp_column rs_id_dbSNP151_GRCh38p7 --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --beta_column beta --output_file ${outpath}/${trait}_${p_type}_${tissue_predix}.GTExv8

done <${LOOPINGFILE}

##############################################################################
##############################################################################

#Read in results and plot


ml R/4.2.0-foss-2020b

cat("SETUP: loading libraries\n")
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

#Add directories to raw path
trait <- 'ANM_BEB'
p_type <- 'exact'
path_spredixcan_here <- '/gpfs/gibbs/pi/huckins/projects/PGCED_CCGWAS/spredixcan'
file_pattern_predixcan <- paste0(trait,'.*',p_type,'.*','GTEXv8')

cat("Reading in S-predixcan output from all tissues  \n")
predixcan_alltissue <- NULL
files_output_spredix <- list.files(path_spredixcan_here, pattern=file_pattern_predixcan)

for (file in files_output_spredix) {
	#file='BMI_2018_Adipose_Subcutaneous.GTEXv8'
	print(file)
	path_here <- paste(path_spredixcan_here,file, sep='/')
	spredixcan_here <- as.data.frame(fread(path_here, header=T, sep=",", quote = "",stringsAsFactors=F))
	#Get tissue from file name, save in df
	tissue_here <- gsub(paste0(trait,'_'), '',file)
	tissue_here <- gsub('.GTEXv8', '',tissue_here)
	tissue_here <- gsub('all\\_', '',tissue_here)
	spredixcan_here$tissue <- tissue_here
	#Add to master
	predixcan_alltissue <- rbind(predixcan_alltissue,spredixcan_here)
}

#Replicate the original significance lines 
#Experiment-wide
p_experiment <- .05/nrow(predixcan_alltissue)
#Tissue
bonf_sig <- predixcan_alltissue %>% group_by(tissue) %>% dplyr::count() %>% mutate(bonf=.05/n) 
p_tissue <- mean(bonf_sig$bonf)

#Left join predixcan with features to get gene start location 
toplot_annot <- left_join(toplot,unique(featuredt[,c('Gene','start','Chr')]), by=c('gene_no_version' = 'Gene'))
#Paste gene start and chr 
toplot_annot$location <- paste(toplot_annot$Chr, toplot_annot$start, sep='.')

#Look at bonferroni to annotate
labeldf <- toplot_annot %>% filter(Is_BMIQTL == 'Overlaps BMI-QTL') %>% filter(pvalue <= p_experiment) %>% mutate(label_here=paste(gene_name,My_Tissue,sep=': '))

#Plot manhattan
pdf(paste(path_output,paste(trait,p_type, 'GTEXv8_Manhattan_exact.pdf', sep='_'), sep='/'), height=10, width=20)
ggplot(toplot_annot,aes(x = as.numeric(location), y = -log10(pvalue), color=as.factor(CHR))) + 
geom_point() + 
theme_classic(base_size=24) +
geom_hline(yintercept=-log10(p_experiment), linetype='dashed', color='red') + 
geom_hline(yintercept=-log10(p_tissue), linetype='dashed', color='darkgray') + 
geom_text_repel(data=labeldf, aes(label=label_here), size=6, max.overlaps=20, color='black') + 
scale_x_continuous(limits=c(1,23),breaks=seq(1.5,22.5), labels=seq(1,22)) + 
labs(x='Chromosome', y='-log10(p)')
dev.off()
