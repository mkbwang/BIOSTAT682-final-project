#!/bin/bash
#SBATCH --job-name=wangmk_jags_logisticregression
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=1g
#SBATCH --cpus-per-task=4
#SBATCH --array=1-11


Rscript binary_classification_jags.R

