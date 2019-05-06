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
if [ "$HOSTNAME" = mpp2-login5 ]; then
    AMOUNT_THREADS=24
    THREAD_PINNING=0-23
else
    AMOUNT_THREADS=28
    THREAD_PINNING=0-27
fi

EXECUTABLE_NAME=inner-loop
export OMP_NUM_THREADS=$AMOUNT_THREADS
mkdir results

likwid-perfctr -g FLOPS_AVX -execpid -C $THREAD_PINNING -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-flops-measurement.out
likwid-perfctr -g FALSE_SHARE -execpid -C $THREAD_PINNING -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-false-sharing-measurement.out
likwid-perfctr -g CACHES -execpid -C $THREAD_PINNING -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-caches-measurement.out
likwid-perfctr -g MEM -execpid -C $THREAD_PINNING -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-memory-measurement.out
likwid-perfctr -g ENERGY -execpid -C $THREAD_PINNING -O -m ./build/$EXECUTABLE_NAME > results/$EXECUTABLE_NAME-energy-measurement.out
