#!/bin/bash

#SBATCH --job-name=PW_run
#SBATCH --output=out.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH -p extended-28core

module load shared
module load R/3.6.2
module load gcc-stack
module load openblas/dynamic/0.2.18
module load lapack/gcc/64/3.6.0
module load JAGS/4.3.0
module load gnu-parallel/6.0

export R_LIBS=/gpfs/home/cyoungflesh/R_libs

cd /gpfs/home/cyoungflesh/penguin_watch_model/Scripts

Rscript 2-model.R

