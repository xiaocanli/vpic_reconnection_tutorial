#!/bin/sh

# The current supported machines include Cori, Frontera, and Perlmutter_cpu or Perlmutter_gpu
MACHINE=$1

if [ "${MACHINE,,}" = "cori" ]; then
    module swap craype-haswell craype-mic-knl
    module unload craype-hugepages2M
    module load cray-hdf5-parallel
elif [ "${MACHINE,,}" = "perlmutter_cpu" ]; then
    module unload craype-hugepages2M
    module load cray-hdf5-parallel
elif [ "${MACHINE,,}" = "perlmutter_gpu" ]; then
    module unload craype-hugepages2M
    module load cray-hdf5-parallel
else
    module load phdf5 # Frontera
fi
module list 2>&1 | tee modules_list

deck_name=reconnection

mkdir spectrum
mkdir particle
mkdir tracer
mkdir field_hdf5
mkdir hydro_hdf5
lfs setstripe -S 8388608 -c 8 spectrum
lfs setstripe -S 8388608 -c 8 particle
lfs setstripe -S 8388608 -c 8 tracer
lfs setstripe -S 8388608 -c 8 field_hdf5
lfs setstripe -S 8388608 -c 8 hydro_hdf5

dir_name=${PWD##*/}
make -f Makefile_${MACHINE,,}
# if [ "${MACHINE,,}" = "cori" ] || [ "${MACHINE,,}" = "perlmutter" ]; then
#     mkdir -p files_sbcast
#     cp $deck_name.Linux files_sbcast
#     cp $(ldd $deck_name.Linux | sed -e "1d" | awk '{print $3}' | paste -s ) files_sbcast
# fi

# # copy the deck to a safe place
# cd ../
# ctime=`date '+%Y%m%d_%H%M%S'`
# tar zcvf ${dir_name}_$ctime.tar.gz $dir_name
# mkdir -p ~/vpic_runs/
# mv ${dir_name}_$ctime.tar.gz ~/vpic_runs/
# cd $dir_name
