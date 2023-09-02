# Modify the Compiling Script
We only provide compiling scripts (cori-knl-openmp-hdf5, perlmutter-openmp-hdf5, tacc-frontera-hdf5) for NERSC's Cori or Perlmutter or TACC's Frontera in this tutorial. However, you can modify the scripts for your platforms without too much effort since the scripts have detailed explanations for each variable—please check out either of the scripts. Usually, you might need to modify the following variables.

* `VCOM`: compiler. It will determine the corresponding compiling flags.
* `VMPI`: MPI implementation. I will determine the modules to load. These modules might have different names on different platforms. Please make corresponding changes.
* `VTHR`: thread model. Either one should work pretty well.
* vector intrinsics: It depends on the CPU used on the platform. Please check the architecture description of your platform to find out which vector intrinsics are supported.
* `set_hdf5`: It is recommended to set it to `ON` for HDF5 IO.

The rest should be similar on different platforms. However, if you find problems, please read the comments in the scripts to debug first.