#!/bin/bash

#PBS -N zeta1
#PBS -A itea_lille-tma4280
#PBS -W group_list=itea_lille-tma4280
#PBS -l walltime=00:01:00
#PBS -l select=2:ncpus=20:mpiprocs=64
#PBS -o output/stdout.txt
#PBS -e output/stderr.txt
#PBS -M asgl@dtu.dk
#PBS -m e

module load GCC/6.4.0-2.28

cd $PBS_O_WORKDIR

mkdir -p -- res
make convtest
