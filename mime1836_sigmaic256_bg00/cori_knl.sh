#!/bin/bash
#
#SBATCH -q regular
#SBATCH -N 64
#SBATCH -t 24:00:00
#SBATCH -C knl,quad,cache
#SBATCH -o vpic%j.out
#SBATCH -e vpic%j.err
#SBATCH -J reconnection
#SBATCH -A m2407
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=phyxiaolee@gmail.com
#SBATCH -L SCRATCH,project

##### These are shell commands
date
module swap craype-haswell craype-mic-knl
module unload craype-hugepages2M
module load cray-hdf5-parallel
# module load lustre-default
module load dws
module list

RUN_DIR=/global/cscratch1/sd/xiaocan/tail_problem/3D_test2

cd $RUN_DIR/files_sbcast

srun mkdir -p /tmp/xiaocanli
for f in ./* ; do sbcast --compress --force --fanout=8 -t 240 $f /tmp/xiaocanli/$f; done
wait
echo "sbcast done!!!"
export LD_LIBRARY_PATH=/tmp/xiaocanli/:$LD_LIBRARY_PATH
ldd /tmp/xiaocanli/reconnection.Linux

###sbcast --compress --force --fanout=8 -t 240 ./untar.sh /tmp/untar.sh &  # Untar after stage in.
#sbcast --compress --force --fanout=8 -t 240 ./reconnection.Linux /tmp/reconnection.Linux &
###sbcast --compress --force --fanout=8 -t 240 ./tar.sh /tmp/tar.sh &
#wait

cd $RUN_DIR

# mkdir spectrum
# mkdir particle
# mkdir tracer
# mkdir field_hdf5
# mkdir hydro_hdf5
# # lfs setstripe -S 16777216 -c 32 field_hdf5
# # lfs setstripe -S 16777216 -c 32 hydro_hdf5
# # lfs setstripe -S 16777216 -c 32 particle
# # lfs setstripe -S 8388608 -c 32 tracer
# # stripe_medium field_hdf5
# # stripe_medium hydro_hdf5
# # stripe_medium particle
# # stripe_medium tracer
# lfs setstripe -S 8388608 -c 24 spectrum
# lfs setstripe -S 8388608 -c 24 particle
# lfs setstripe -S 8388608 -c 24 tracer
# lfs setstripe -S 8388608 -c 24 field_hdf5
# lfs setstripe -S 8388608 -c 24 hydro_hdf5

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#Quad Cache:
export PMI_MMAP_SYNC_WAIT_TIME=400
time srun -n 4096 -N 64 -c 4 --cpu_bind=cores /tmp/xiaocanli/reconnection.Linux --tpp 4

date
echo 'Done'
