#!/bin/bash

#SBATCH --job-name=boot2d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=5300
#SBATCH --time=24:00:00
#SBATCH --output=boot2d.out



# --ntasks-per-node will be used in doParallel.R to specify the number of 
# cores to use on the machine. Using 16 will allow us to use all cores
# on a sandyb node

module unload R
module load netcdf/4.4.1 R/3.3
cd /project/moyer/mahaugen/r_code

Rscript boot2d.R