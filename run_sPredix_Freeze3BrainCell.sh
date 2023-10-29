#!/bin/bash

# load environment
ml miniconda
conda activate py3_envNEW_spredix

# make directory for project
mkdir /home/ekw28/AN-CellTypeS-Predixcan/AN_Freeze3_brainCellTypes
mkdir /home/ekw28/AN-CellTypeS-Predixcan/AN_Freeze3_brainCellTypes/scratch                                                            



## set software file locations
SPREDIX_path=/gpfs/gibbs/pi/huckins/software/MetaXcan/software/SPrediXcan.py
outpath=/home/ekw28/AN-CellTypeS-Predixcan/AN_Freeze3_brainCellTypes
scratch=/home/ekw28/AN-CellTypeS-Predixcan/AN_Freeze3_brainCellTypes/scratch             
cd $scratch

## change this to suit GWAS
trait=AN
trait_filename=freeze3.AN.GWAS.sumstats.clean
path_GWAS=/home/ekw28/AN-CellTypeS-Predixcan

## set locations of predictor models
# first cell type
cellType_predix= #place tissue here

PREDICT_DB_path=/gpfs/gibbs/pi/huckins/girgenti_celltype_models/${cellType_predix}.db
COVARIANCE_PATH=/gpfs/gibbs/pi/girgenti_celltype_models/${cellType_predix}.txt.gz


##Change the snp_column, effect_allele_column, non_effect_allele_column, pvalue_column, or_column to match the colnames in the GWAS
#https://github.com/hakyimlab/MetaXcan/wiki/Command-Line-Reference#spredixcanpy
python ${SPREDIX_path} --model_db_path ${PREDICT_DB_path} --covariance ${COVARIANCE_PATH} --gwas_folder ${path_GWAS} --gwas_file_pattern ${trait_filename} --snp_column rs.id.dbSNP151.GRCh38p7 --effect_allele_column effect_allele --non_effect_allele_column noneffect_allele --pvalue_column p --or_column OR --output_file ${outpath}/${trait}_${p_type}_${cellType_predix}.carina_cell_type

