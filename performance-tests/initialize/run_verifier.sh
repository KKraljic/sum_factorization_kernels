#!/bin/bash

#SBATCH -o ./flops-out.txt
#SBATCH -D .

#SBATCH -J LIKWID-K
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --export=NONE
#SBATCH --time=00:05:00

# LOAD MODULE
module load likwid/4.3
#AMOUNT_THREADS=28
#export OMP_NUM_THREADS=$AMOUNT_THREADS
mkdir results
likwid-perfctr -g FLOPS_AVX -execpid -C 0-27 -O -m build/verifier > results/icc_flops.out
