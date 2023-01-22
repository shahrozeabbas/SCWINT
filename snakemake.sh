#!/bin/sh

#SBATCH --mem=600g
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=36
#SBATCH --partition=largemem
#SBATCH --mail-type=ALL,TIME_LIMIT_50

source /data/abbass2/Apps/conda/bin/activate snakes

snakemake --use-conda --cores 32