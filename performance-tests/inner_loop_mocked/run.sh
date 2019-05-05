#!/bin/bash

#SBATCH -o ./inner-loop-out.txt
#SBATCH -D .

#SBATCH -J LIKWID-K
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --export=NONE
#SBATCH --time=00:60:00

# LOAD MODULE
module load likwid/4.3
AMOUNT_THREADS=28
EXECUTABLE_NAME=inner-loop
export OMP_NUM_THREADS=$AMOUNT_THREADS
mkdir results

likwid-perfctr -g FLOPS_AVX -execpid -C 0-27 -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-flops-measurement.out
likwid-perfctr -g FALSE_SHARE -execpid -C 0-27 -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-false-sharing-measurement.out
likwid-perfctr -g CACHES -execpid -C 0-27 -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-caches-measurement.out
likwid-perfctr -g MEM -execpid -C 0-27 -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-memory-measurement.out
likwid-perfctr -g ENERGY -execpid -C 0-27 -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-energy-measurement.out
