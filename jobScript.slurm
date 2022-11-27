#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --account=PCS0268
#SBATCH --job-name=AL_FINAL

N=50
export OMP_NUM_THREADS=20
cd $SLURM_SUBMIT_DIR/Default
make all

# test run with debug
./Household.exe  30 10 1 > /dev/null


echo "N p time" > results.csv
for p in 1, 2, 3, 4, 5, 10, 20
do
    ./Household.exe ${N} ${p} >> results.csv
done

