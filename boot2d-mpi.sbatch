#!/bin/bash


#SBATCH --ntasks=240
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=6300
#SBATCH --output=boot2d-mpi.out
#SBATCH --time=8:05:00

module unload R
module load netcdf/4.4.1 R/3.3
cd /project/moyer/mahaugen/r_code

mpirun -np 240 Rscript boot2d-mpi.R