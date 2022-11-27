#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --account=PCS0268
#SBATCH --job-name=AL_FINAL

export OMP_NUM_THREADS=20
cd $SLURM_SUBMIT_DIR/Default
make all

# test run with debug
./Household.exe  30 10 1 > /dev/null

N=50

echo "N p tdTime evTime" > results.csv
for N in 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
do
    for p in 1, 2, 3, 4, 5, 10, 20
    do
        ./Household.exe ${N} ${p} >> results.csv
    done
done

