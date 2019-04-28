#!/bin/bash

#SBATCH -o ./matvec-job.out
#SBATCH -D .

#SBATCH -J LIKWID-K
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --export=NONE
#SBATCH --time=00:30:00

# LOAD MODULE
module load likwid/4.3
mkdir results

likwid-perfctr -g FLOPS_AVX -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-flops-measurement.out
likwid-perfctr -g FALSE_SHARE -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-false-sharing-measurement.out
likwid-perfctr -g CACHES -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-caches-measurement.out
likwid-perfctr -g MEM -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-memory-measurement.out
likwid-perfctr -g ENERGY -execpid -C 1 -O -m ./build/matvec-measurement > results/matvec-energy-measurement.out




