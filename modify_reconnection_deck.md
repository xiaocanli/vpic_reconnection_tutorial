# Modify reconnection deck
If you decide to modify the reconnection deck (`reconnection.cc`), you need to make sure that the parameters are reasonable (e.g., the box size, particle number, and memory usage). There is a `run_checklist.txt` for checking the deck before running the simulation. The source code includes detailed comments on all the variables. Some of them might be outdated, so please help to improve the descriptions. Some variables that we usually change are listed below.
  - `particle_tracing`: 0 or 1
  - `particle_select`: make sure that the total number of tracked particles is not too large (about 1 million)
  - `mi_me`: mass ratio
  - `vthe`: thermal speed.
  - `bg`: guide field strength. $10^{-6}$ indicates 0.0. `0.1` indicates 0.1.
  - `taui`: simulation time in $\Omega_{ci}^{-1}$.
  - `nppc`: number of particles/cell/species. 100 is usually enough.
  - `nx`, `ny`, `nz` are the simulation grid numbers along each dimension. `ny=1` for 2D simulations.
  - `Lx`, `Ly`, `Lz` are simulation box sizes along each dimension.
  - `topology_x`, `topology_y`, `topology_z`: MPI topology.
  - `interval`: reasonable number to generate a few hundreds of frames.
  - `tracer_interval`: reasonable number so that the number of tracer time steps is about 1000.
  - `eparticle_interval`: dump all particles. A few times of `interval`
  - `nx_zone`, `ny_zone`, `nz_zone`: see the description in the deck.
  - `stride_particle_dump`: dump one particle every a few particles

  To dump particle data, you need to set `eparticle_interval` and `Hparticle_interval` to a smaller value (the same as or maybe a few times of the fields_interval). You also need to uncomment the lines
  ```cpp
  sprintf(subdir,"particle/T.%d/eparticle",step());
  dump_particles("electron",subdir);
  sprintf(subdir,"particle/T.%d/hparticle",step());
  dump_particles("ion",subdir);
  ```
  Also, make sure to change `stride_particle_dump` if you only need a fraction of all the particles. For example, you can set `stride_particle_dump = 20` to dump 5\% of all the particles. In this method, particle data will be in binary format, and each MPI rank will write one particle file/step. This kind of N->N IO is not allowed in some HPC platforms (e.g., Frontera). That's why you may need to write the particle data in HDF5 format. Checkout `dump_with_h5part.cc` on the procedure. If you decide to do that, we can discuss the details.