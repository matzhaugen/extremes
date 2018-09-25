#!/bin/bash


#SBATCH --ntasks=600
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=6300
#SBATCH --output=boot2d-mpi.out
#SBATCH --time=20:05:00

module unload R
module load netcdf/4.4.1 R/3.3
cd /project/moyer/mahaugen/r_code

mpirun -np 600 Rscript mainBig-mpi.R