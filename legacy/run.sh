#!/bin/bash

#SBATCH -o ./dg-integrate.txt
#SBATCH -D .

#SBATCH -J LIKWID-DG-INT-K
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --export=NONE
#SBATCH --time=00:00:30

# LOAD MODULE
module load mpi.intel
module load likwid/4.3
AMOUNT_THREADS=28
export OMP_NUM_THREADS=$AMOUNT_THREADS

POLYNOMIAL_DEGREE=4
PROBLEM_SIZE=10000000
N_ITERATIONS=200
MODE=1

RESULTS=results
mkdir $RESULTS
likwid-perfctr -g CACHES -execpid -C 0-27 -O -m build/test_dg_integrate $POLYNOMIAL_DEGREE $PROBLEM_SIZE $N_ITERATIONS $MODE > $RESULTS/dg-integrate-caches.out
likwid-perfctr -g FLOPS -execpid -C 0-27 -O -m build/test_dg_integrate $POLYNOMIAL_DEGREE $PROBLEM_SIZE $N_ITERATIONS $MODE > $RESULTS/dg-integrate-flops.out
