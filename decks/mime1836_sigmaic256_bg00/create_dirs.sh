#!/bin/sh

rm -rf spectrum
rm -rf particle
rm -rf tracer
rm -rf field_hdf5
rm -rf hydro_hdf5

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
