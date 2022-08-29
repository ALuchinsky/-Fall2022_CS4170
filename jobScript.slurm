#!/bin/bash

#SBATCH --ntasks=2 
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --account=PCS0253
#SBATCH --job-name=JOB_NAME

export OMP_NUM_THREADS=12
cd $SLURM_SUBMIT_DIR/Default
make all
./OpenMP > results.csv
