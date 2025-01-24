#!/bin/sh

# MACHINE=Frontera
MACHINE=Perlmutter

if [ "$MACHINE" = "Perlmutter" ]; then
    module load cpu
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
mkdir fields-avg-hdf5
mkdir hydro-avg-hdf5
lfs setstripe -S 33554432 -c 16 spectrum
lfs setstripe -S 33554432 -c 16 particle
lfs setstripe -S 33554432 -c 16 tracer
lfs setstripe -S 33554432 -c 16 field_hdf5
lfs setstripe -S 33554432 -c 16 hydro_hdf5
lfs setstripe -S 33554432 -c 16 fields-avg-hdf5
lfs setstripe -S 33554432 -c 16 hydro-avg-hdf5

dir_name=${PWD##*/}
if [ "$MACHINE" = "Perlmutter" ]; then
    make -f Makefile_perlmutter
else
    make -f Makefile_tacc # Frontera
fi

# # copy the deck to a safe place
# cd ../
# ctime=`date '+%Y%m%d_%H%M%S'`
# tar zcvf ${dir_name}_$ctime.tar.gz $dir_name
# mkdir -p ~/vpic_runs/
# mv ${dir_name}_$ctime.tar.gz ~/vpic_runs/
# cd $dir_name
