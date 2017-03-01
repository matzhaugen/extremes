#!/bin/bash

#SBATCH --job-name=main2d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --partition=bigmem2
#SBATCH --time=04:00:00



# --ntasks-per-node will be used in doParallel.R to specify the number of 
# cores to use on the machine. Using 16 will allow us to use all cores
# on a sandyb node

module load netcdf/4.4.1 R/3.2 openmpi/2.0.1   
cd /project/moyer/mahaugen/r_code

Rscript main2d.R