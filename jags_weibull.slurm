#!/bin/bash
#SBATCH --job-name=wangmk_jags_weibull
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=1g
#SBATCH --cpus-per-task=4
#SBATCH --array=1-11


Rscript jags_weibull.R
