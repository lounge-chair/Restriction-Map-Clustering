#!/bin/bash
#SBATCH --qos=boucher
##SBATCH --partition=hpg2-compute
#SBATCH --job-name=cluster
#SBATCH --mail-user=kingdgp.lfc@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output output
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --time=4-00:00:00

/usr/bin/time -v Rscript clustering.r