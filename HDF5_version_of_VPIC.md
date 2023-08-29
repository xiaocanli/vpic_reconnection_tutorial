# HDF5 version of VPIC

Check out the hdf5_rebase branch of VPIC

```bash
git clone https://github.com/lanl/vpic
cd vpic
git checkout hdf5_rebase
```
> **_NOTE:_**  The hydro dump in the `hdf5_rebase` branch <s>https://github.com/lanl/vpic</s> has problems (https://github.com/lanl/vpic/issues/146). Please use the `hdf5_rebase` branch in https://github.com/xiaocanli/vpic instead.

I usually rename the directory as `vpic_hdf5_rebase` to differentiate with other branches.

Copy `cori-knl-openmp-hdf5` to `vpic_hdf5_rebase/arch`

```bash
cd vpic_hdf5_rebase
mkdir build
cd build
../arch/cori-knl-openmp-hdf5
```

On Perlmutter, for the CPU nodes, use `perlmutter-openmp-hdf5` instead. For the GPU nodes, please refer [decks/problems.md](decks/problems.md) for an instruction.

On Frontera, use `../arch/tacc-frontera-hdf5` instead. It will produce the `bin/vpic` for compiling the deck. For Frontera, you need to modify `bin/vpic` to link the HDF5 library by adding `Wl,-rpath,$TACC_HDF5_LIB -L$TACC_HDF5_LIB -lhdf5 -lz` to `bin/vpic`. An example is given in this directory.

Before compiling the deck, you need to load `module load cray-hdf5-parallel` on Cori or `module load phdf5` on Frontera. Modify `VPIC_DIR` in `Makefile_cori` or `Makefile_tacc`  and compile the deck using `make -f Makefile_cori` or `make -f Makefile_tacc`. There is a script `compile_deck.sh` for helping compile the deck and copy the source files to a safe place. 

To do a test on an interactive session on Cori,

```bash
salloc -N 1 -q interactive -C knl,quad,cache -t 04:00:00 -L SCRATCH
./run_test.sh
```

On Frontera, you can request an interactive session using

```bash
idev -p normal -N 2 -n 64 -m 60 -A PHY20020
ibrun ./reconnection.Linux
```
The code will output fields and hydro data to `field_hdf5` and `hydro_hdf5`, respectively. Note that electron hydro data and ion hydro data are saved in different files. You can check what is included in the `.h5` files using `h5dump` after loading the HDF5 modules. To learn more about these command-line tools, please check out [HDF5 TUTORIAL: COMMAND-LINE TOOLS FOR VIEWING HDF5 FILES](https://support.hdfgroup.org/HDF5/Tutor/cmdtoolview.html). The whole [HDF5 TUTORIAL](https://support.hdfgroup.org/HDF5/Tutor/) is very helpful if you are not familiar with HDF5. For example, `h5dump -H fields_0.h5` gives
```
HDF5 "fields_0.h5" {
GROUP "/" {
   GROUP "Timestep_0" {
      DATASET "cbx" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
      DATASET "cby" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
      DATASET "cbz" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
      DATASET "ex" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
      DATASET "ey" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
      DATASET "ez" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 256, 1, 256 ) / ( 256, 1, 256 ) }
      }
   }
}
}
```
`cbx`, `cby`, `cbz` are magnetic fields and `ex`, `ey`, `ez` are electric fields. Similarly, the hydro files include `jx`, `jy`, and `jz` for current density, `ke` for kinetic energy density, `px`, `py`, and `pz` for momentum density, `rho` for charge density, `txx`, `txy`, `tyy`, `tyz`, `tzx`, and `tzz` for stress tensor. These are similar to the data before `translate`.
- `jx/rho` will be `vex` or `vix`
- `|rho|` will be `ne` or `ni`
- `px/|rho|/particle_mass` with be the commonly used `uex` or `uix`. `particle_mass` is 1 for electrons or `mi_me` for ions.
- `txx - (jx/rho)*px` will be `pexx` or `pixx`. `txy - (jx/rho)*py` will be `pexy` or `pixy`. `txy - (jy/rho)*px` will be `peyx` or `piyx`. Similar for other pressure tensor components.

You can use [quick_check_vpic](https://github.com/xiaocanli/quick_check_vpic) to check the data. For detailed analysis using the HDF5 files, please follow `hdf5_analysis_example.py` in this directory. It gets the VPIC simulation information and plot $j_y$. Reading other variables is similar. After reading into the memory, the rest will be the same. You can still use most of your analysis code. If the data is really large (e.g., in 3D simulations), we use ParaView or VisIt to visualize the data. In `field_hdf5/` and `hydro_hdf5/`, there are files ending with `.xdmf`, which can be loaded into ParaView or VisIt directly. There is a python script `mime1836_sigmaic256_bg00/gen_xdmf.py` to generate a single `xdmf` file for both fields and hydro data when necessary.

For 2D simulations with periodic boundary conditions along $x$ and conducting boundaries along $z$, you can use `ay_gda_hdf5.f90` in the deck to calculate the $y$-component of the vector potential $A_y$ from the magnetic field data. Then, you can plot the magnetic field lines the contour of $A_y$. Please change a few parameters before compiling the code: `it1`, `it2`, `interval`, `hdf5_fields`, `nx`, `nz`, `xmax`, and `zmax`. These parameters can be obtained from `info` and file information in `field_hdf5/`. On Cori, you can compile the code using
```sh
module load cray-hdf5 cray-fftw
ftn -o ay_gda_hdf5 ay_gda_hdf5.f90
```
Then, you can run the code `./ay_gda_hdf5`, which will generate a file `data/Ay.gda`. The file is in binary format and has all time frames in the same file (`nt = it2-it1+1` frames). The data is organized in memory `nx*nz*nt`, where `nx` changes the fastest.
