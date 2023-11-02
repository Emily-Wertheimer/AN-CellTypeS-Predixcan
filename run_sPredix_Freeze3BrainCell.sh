######################################################################################
## submit job script text
#!/bin/bash
#SBATCH -J --job-name=run_sPredix_excitatoryNeurons_AN
#SBATCH --mem=20G
#SBATCH --cpus-per-task=40
#SBATCH --partition=day
#SBATCH --time=1-

## load modules
#module load #module here

#script name
runSpredix_freeze3AN_brainCell.sh

######################################################################################
## location of software & environments
SPREDIX_path=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py
module purge && module load miniconda
#source=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/conda_env.yaml
conda activate py3_env_spredix

# Set software file locations
outpath=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo
scratch=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/scratch
cd /vast/palmer/home.mccleary/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo

# Set GWAS trait information
trait=AN
trait_filename=freeze3.AN.GWAS.sumstats.clean
path_GWAS=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo #/freeze3.AN.GWAS.sumstats.clean

######################################################################################
### submit 1 tissue per job
tissue_predix=BasoAmyg_Excitatory_Neurons # options are: BasoAmyg_, dACC_, DLPFC_, MedialAmyg_
PREDICT_DB_path=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/Excitatory_Neurons/BasoAmyg_Excitatory_Neurons_weights.db
COVARIANCE_PATH=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/BasoAmyg_Excitatory_Neurons_LD_covariance_matrix.txt.gz

## Specify GWAS column names (modify these based on your GWAS data)
snp_column="rs.id.dbSNP151.GRCh38p7"
effect_allele_column="effect_allele"
non_effect_allele_column="noneffect_allele"
pvalue_column="p"
or_column="OR"

## Run spredixcan
#python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename} --snp_column rs.id.dbSNP151.GRCh38p7 --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --or_column OR --output_file ${outpath}/${trait}_${tissue_predix}.carina.cell.type​
python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename} --snp_column ${snp_column} --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --or_column OR --output_file ${outpath}/${trait}_${tissue_predix}.carina.cell.type​

######################################################################################

### loop thru tissues in one job (can't get this to work") 
## Set locations of predictor models (provide an actual value for cellType_predix)
#tissueType_predix=("BasoAmyg_" "dACC_" "DLPFC_" "MedialAmyg_")
#cellType_predix=Excitatory_Neurons 

#PREDICT_DB_path=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo #/${tissue}${cellType_predix}.db
#COVARIANCE_PATH=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo #/${tissue}${cellType_predix}.txt.gz

## Specify GWAS column names (modify these based on your GWAS data)
#snp_column="rs.id.dbSNP151.GRCh38p7"
#effect_allele_column="effect_allele"
#non_effect_allele_column="noneffect_allele"
#pvalue_column="p"
#or_column="OR"

## Run SPrediXcan for e/ tissue type
#LOOPINGFILE=/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/Excitatory_Neurons
#cd $scratch
#for tissue in ${tissueType_predix[@]}; do
 #   python ${SPREDIX_path} \
 #   --model_db_path ${PREDICT_DB_path}/${tissue}/${cellType_predix}.db \
  #  --covariance ${COVARIANCE_PATH} ${tissue} ${cellType_predix}.txt.gz \
   # --gwas_folder ${path_GWAS} \
    #--gwas_file_pattern ${trait_filename} \
   # --snp_column ${snp_column} \
   # --effect_allele_column ${effect_allele_column} \
   # --non_effect_allele_column ${non_effect_allele_column} \
   # --pvalue_column ${pvalue_column} \
   # --or_column ${or_column} \
   # --output_file ${outpath}/${trait}_${tissueType_predix}${cellType_predix}.carina_cell_type
#done <${LOOPINGFILE}
