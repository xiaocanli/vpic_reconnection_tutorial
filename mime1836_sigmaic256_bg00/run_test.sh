#!/bin/bash
date
module swap craype-haswell craype-mic-knl
module unload craype-hugepages2M
module load cray-hdf5-parallel
module list

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#Quad Cache:
# time srun -n 64 -N 1 -c 4 --cpu_bind=cores ./reconnection.Linux --tpp 4
time srun -n 64 -N 1 -c 4 --cpu_bind=cores ./reconnection.Linux --tpp 4 --restore restore1/restore.0

date
echo 'Done'


