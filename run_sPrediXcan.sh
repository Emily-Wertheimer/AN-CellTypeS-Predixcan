####################################################################
# Original author: Rebecca Signer, MS 
# Modified by: Emily Wertheimer 

# Hakyim lab MetaXcan github repo: https://github.com/hakyimlab/MetaXcan/tree/master​

####################################################################
# get connected to hpc 
ssh ekw28@mccleary.ycrc.yale.edu

#Request more memory/time
salloc -t 6:00:00 --mem=32G #tailor to what you need

####################################################################
## load conda env to run spredix 

ml miniconda

# first time making environment
module purge && module load miniconda
conda env create -n py3_env_spredix --file /gpfs/gibbs/pi/huckins/software/MetaXcan/software/conda_env.yaml
module purge && module load miniconda
conda activate py3_env_spredix

# every time after
#conda activate py3_env_spredix # env loaded w/ .yaml 
conda activate py3_envNEW_spredix # manually loaded env

####################################################################
#Location of software
SPREDIX_path=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py

​#Define outpath and scratch
outpath=/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/data
scratch=/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/scratch #<SAVE .e/.o file LOCATION>
cd $scratch
​
#Change this to suite your GWAS
trait=AN_GWAS_2019
trait_filename=anGWAS2019_sumstats_clean_2.txt

Define *clean* GWAS path and columns 
path_GWAS=/gpfs/gibbs/project/levy_ifat/ekw28/sPrediXcan_AN/data

####################################################################
#Define location of gtex v8 predictor model 

# first tissue 
tissue_predix=Adipose_Subcutaneous

PREDICT_DB_path=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_Adipose_Subcutaneous.db
#PREDICT_DB_path=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_${tissue_predix}.db
COVARIANCE_PATH=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_Adipose_Subcutaneous.txt.gz
#COVARIANCE_PATH=/gpfs/gibbs/pi/huckins/resources/predixcan/GTEX_v8/elastic_net_models/en_{tissue_predix}.txt.gz

​
#https://github.com/hakyimlab/MetaXcan/wiki/Command-Line-Reference#spredixcanpy
#Change the snp_column, effect_allele_column, non_effect_allele_column, pvalue_column, or_column to match the colnames in the GWAS
​
python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename} --snp_column rs.id.dbSNP151.GRCh38p7 --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --or_column OR --output_file ${outpath}/${trait}_${p_type}_${tissue_predix}.GTExv8​
