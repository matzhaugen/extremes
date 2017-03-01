#!/bin/sh

#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=5200
#SBATCH --output=hello.out
#SBATCH --time=01:05:00

mpirun -np 24 Rscript hello.R