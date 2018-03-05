#!/bin/bash

#PBS -N pi
#PBS -A itea_lille-tma4280
#PBS -W group_list=itea_lille-tma4280
#PBS -l walltime=00:01:00
#PBS -l select=2:ncpus=20:mpiprocs=64
#PBS -o output/stdout.txt
#PBS -e output/stderr.txt
#PBS -M asgl@dtu.dk
#PBS -m e

cd $PBS_O_WORKDIR
module load gcc openmpi

mkdir -p -- res
make convtest
