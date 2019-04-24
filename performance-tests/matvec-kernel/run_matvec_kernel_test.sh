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
mkdir results

#likwid-perfctr -g FLOPS_AVX -execpid -C 0-27 -O -m build/matvec-measurement > results/matvec-measurement.out
likwid-perfctr -g FLOPS_AVX -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-measurement.out

