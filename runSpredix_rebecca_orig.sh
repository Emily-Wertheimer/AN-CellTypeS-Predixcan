########################################
#Location of software
SPREDIX_path=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py
ml python
####################################################################
​
#Define outpath and scratch
outpath= #<SAVE LOCATION>
scratch= #<SAVE .e/.o file LOCATION>
cd $scratch
​
#Change this to suite your GWAS
traits=("IBD_2017" "Asthma_2022" "MS_2019" "T1D_2021")
trait_filenames=("IBD_2017-cd_build37_40266_20161107-OR_StdErr.txt.gz" "Asthma_2022-Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1-OR_se.txt.gz" "MS_2019-discovery_metav3.0-OR_se.txt.gz" "T1D_2021-GCST90014023_buildGRCh38-OR_standard_error.txt.gz")
​
#Define GWAS path and columns 
#Raw AN phase 2 is /gpfs/gibbs/pi/huckins/resources/gwas_spredixcan_sumstats/PGC/an2019/gwas_raw
#You need to clean the GWAS so you know the build matches the CMC model build and the effect allele and stat columns are determined
path_GWAS= #<Path to your cleaned GWAS>
​
####################################################################
#Define location of CMC predictor model (can be downloaded from online) - check with Laura, I only use GTEx which is located
tissue_predix=DLPFC
PREDICT_DB_path=/gpfs/gibbs/pi/huckins/resources/predixcan/CMC/CMC_DLPFC-PrediXcan_Models/${tissue_predix}_newMetax.db
COVARIANCE_PATH=/gpfs/gibbs/pi/huckins/resources/predixcan/CMC/CMC_DLPFC-PrediXcan_Models/${tissue_predix}.cov.txt.gz
jobname=${trait_here}_${tissue_predix}_spredixcan
​
#https://github.com/hakyimlab/MetaXcan/wiki/Command-Line-Reference#spredixcanpy
#Change the snp_column, effect_allele_column, non_effect_allele_column, pvalue_column, or_column to match the colnames in the GWAS
​
bsub -q express -P acc_psychgen -J ${jobname} -n 1 -W 1:00 -R rusage[mem=20000] -o ${jobname}.o -e ${jobname}.e python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename_here} --snp_column rs_id_dbSNP151_GRCh38p7 --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --or_column OR --output_file ${outpath}/${trait_here}_${tissue_predix}.CMC
​
​
done <${LOOPINGFILE}
