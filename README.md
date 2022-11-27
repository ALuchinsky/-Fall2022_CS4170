# CS5170, Final Project
Aleksei Luchinsky

## Installation

This project was created from Dr. Green's OpenMP template. To build to sources you should simply do

    cd Default
    make all
from the the root directory of the repository. As a result, the following executable will be created:

* serial_Jacoby.exe:        Serial version of the Jacoby transformation algorithm, works
* parallel_Jacoby.exe       Parallel version of the Jacoby transformation algorithm, does not work
* Power.exe                 An attempt to use the Power algorithm, does not work
* Household.exe             Three-Diagonalization using household transformations

Only the last one will be used in the following, the others are kept just for history.

The `Household.exe` program can be run as follows:

    ./Household.exe <N=10> <p=10> <debug=0>
where `N` is the number of rows/columns is the matrix, `p` is the desired number of threads, and `debug` is a flag, that specifies wither to save some debug files.

The program will create a random $N \times N$ matrix, convert it into 3-diagonal form using a set of Household transformations, and then search for eigenvalues using `Bisect` method. As a result, the following line will be printed to the standard output:

    <N> <p> <tdTime> <evTime>
Here `N` and `p` are matrix dimension and number of threads, while `tdTime` and `evTime` are times, required to make Household transformation and calculate the eigenvalues respectively.
If the option `debug=1` was used, initial random matrix, its 3-diagonalized version and a set of eigenvalues will be saved to the files  `in_hA.txt`,   `out_hA.txt`, and `ev_out.txt`.

If you want to do to the preforance analysis, you can run the program on OSC cluster in batch mode. To do this you should go back to the repository's home directory and run the command

     sbatch jobScript.slurm
This script will run the main program with several sets of `N`, `p` configurartions, and save the statistics to `Default/results.csv` file. This file was alayzed by the `analysis/final_analysis.ipynb` notebook.

## Analysis




