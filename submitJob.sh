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

