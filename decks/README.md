# Problem Decks
## Simulation decks for the CPU nodes
* [2D-forcefree](2D-forcefree): nonrelativistic reconnection starting from a forcefree current sheet
* [mime1836_sigmaic256_bg00](mime1836_sigmaic256_bg00): relativistic reconnection
* [2D-Lx150-bg0.2-150ppc-beta002](2D-Lx150-bg0.2-150ppc-beta002): 2D nonrelativistic reconnection
* [3D-Lx150-bg0.2-150ppc-beta002](3D-Lx150-bg0.2-150ppc-beta002): 3D nonrelativistic reconnection. It is similar to the 2D deck.
* [mime1_sigmae40_vthe05_db1](mime1_sigmae40_vthe05_db1): relativistic turbulence
* [sigmae100_bg00_400de](sigmae100_bg00_400de): relativistic reconnection with tracer particles. Please follow `tracer_eb.ipynb` on how to run the simulation and analyze the data.

## Reconnection using vpic-kokkos for GPU nodes

> [!WARNING]
> This part is outdated.

VPIC can be run on GPUs using the [vpic-kokkos](https://github.com/lanl/vpic-kokkos). The official documentation is at [VPIC-Kokkos’s documentation](https://lanl.github.io/vpic-kokkos/index.html). Please read it and maybe do some small tests. We are still working on migrating the simulations decks to vpic-kokkos, so the following instruction might change. We added HDF5 IO for vpic-kokkos and included tracer particles, which require some modifications of the devel branch of vpic-kokkos. Please use the modified vpic-kokkos at https://github.com/xiaocanli/vpic-kokkos
```sh
git clone --recursive https://github.com/xiaocanli/vpic-kokkos
cd vpic-kokkos
mkdir build
cd build
../arch/Perlmutter_GPU
```
which will generate `bin/vpic` in the `build/` directory. Since tracer particles are included in the simulation deck, we need to modify `bin/vpic` for compiling the deck with tracers. Copy [reconnection_kokkos](reconnection_kokkos) to `$SCRATCH` for running the simulations. Then, copy `bin/vpic` to `reconnection_kokkos` and open `reconnection_kokkos/vpic` with the editor of your choice.
* Find `KOKKOS_CONTAINER_INCLUDES` and add
    ```sh
    KOKKOS_ALGORITHM_INCLUDES=-I/global/u2/x/xiaocan/vpic_sources/vpic-kokkos/kokkos/algorithms/src
    echo $KOKKOS_CONTAINER_INCLUDES
    ```
    below it. **Make sure to change `/global/u2/x/xiaocan/vpic_sources/vpic-kokkos` to where your `vpic-kokkos` is.** Then, add `$KOKKOS_ALGORITHM_INCLUDES` between `$KOKKOS_CONTAINER_INCLUDES` and `$KOKKOS_LIBS` in the bottom two long lines.
* Add `-DVPIC_ENABLE_HDF5` behind `bin/CC` and `bin/nvcc_wrapper` in the bottom two long lines.
* Add `test_particle/advance_p_test_particle.cc test_particle/boundary_p_test_particle.cc test_particle/move_p_test_particle.cc` after `-DINPUT_DECK='"'$1'"'` in the bottom two long lines.

The deck `reconnection_kokkos/` includes a file `vpic_modified` for reference. With the modified `vpic`, you can then compile the reconnection deck using
```sh
./compile_deck.sh perlmutter_gpu
```
which will generate an executable `reconnection.Linux`. You can then request an interactive node to run a test, for example,
```sh
salloc --nodes 1 --qos interactive --time 04:00:00 --constraint gpu --account=m4054_g
module load gpu cray-hdf5-parallel
srun -n 4 -G 4 ./reconnection.Linux --tpp 4 –kokkos-num-devices=4
```
which uses 4 MPI ranks and 4 GPUs on the node to run the simulation. It will take a few minutes. The simulation generates the same files (`field_hdf5`, `hydro_hdf5`, `spectrum`, `tracer`) as the previous CPU runs.
