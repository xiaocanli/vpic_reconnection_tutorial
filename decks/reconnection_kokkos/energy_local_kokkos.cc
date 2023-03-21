{
  // energy band diagnostics
  static int energy_diagnostics_first_time = 1;
  static std::vector<edata *> edParams; // vector to hold species parameters, such as vth

  if (should_dump(spectrum)) {
    if (energy_diagnostics_first_time) {
      sim_log("initializing the energy diagnostics");

      edParams.push_back(&global->ede);
      edParams.push_back(&global->edH);

      energy_diagnostics_first_time = 0;
    } // first time called

    char fname[256];
    char spectrum_scratch[128];
    char subspectrum_scratch[128];
    char dset_name[64];
    sprintf(spectrum_scratch, "spectrum");
    dump_mkdir(spectrum_scratch);
    sprintf(subspectrum_scratch, "%s/T.%d/", spectrum_scratch, step());
    dump_mkdir(subspectrum_scratch);

    int nsp = global->nsp;

    // Prepare for HDF5 IO
    sprintf(fname, "%s/spectrum_%d.h5", subspectrum_scratch, step());

    double el1 = uptime();

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(fname , H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    int nbins = global->nbins;
    float emin_spect_log = log10(global->emin_spect);
    float dloge = (log10(global->emax_spect) - emin_spect_log) / (nbins - 1);
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    int nx2 = nx + 2;
    int ny2 = ny + 2;
    int nxy2 = nx2 * ny2;
    int nx_zone = global->nx_zone;
    int ny_zone = global->ny_zone;
    int nz_zone = global->nz_zone;
    int nzones_x = (grid->nx + nx_zone - 1) / nx_zone;
    int nzones_y = (grid->ny + ny_zone - 1) / ny_zone;
    int nzones_z = (grid->nz + nz_zone - 1) / nz_zone;
    float incells = 1.0 / static_cast<float>(global->nx_zone * global->ny_zone * global->nz_zone);

    // Prepare to write the spectrum data 
    int ix, iy, iz ;
    RANK_TO_INDEX( int(rank()), ix, iy, iz );
    hsize_t global_sizes[4], local_sizes[4], offsets[4];

    global_sizes[0] = global->topology_z * nzones_z;
    global_sizes[1] = global->topology_y * nzones_y;
    global_sizes[2] = global->topology_x * nzones_x;
    global_sizes[3] = nbins + 3;  // including 3 components of B-field
    local_sizes[0] = nzones_z;
    local_sizes[1] = nzones_y;
    local_sizes[2] = nzones_x;
    local_sizes[3] = nbins + 3;
    offsets[0] = iz * nzones_z;
    offsets[1] = iy * nzones_y;
    offsets[2] = ix * nzones_x;
    offsets[3] = 0;

    hid_t filespace = H5Screate_simple(4, global_sizes, NULL);
    hid_t memspace =  H5Screate_simple(4, local_sizes, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);

    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, local_sizes, NULL);
    el1 = uptime() - el1;
    sim_log("TimeHDF5Open: "<< el1 << " s");  //Easy to handle results for scripts

    // Average magnetic field in local zones
    typedef Kokkos::View<float***[3]>  View3DType;
    View3DType bfield_avg_d( "magnetic_field_avg", nzones_z, nzones_y, nzones_x);
    View3DType::HostMirror bfield_avg_h = Kokkos::create_mirror_view( bfield_avg_d );
    Kokkos::deep_copy(bfield_avg_d, 0.0f);
    auto& k_interp = interpolator_array->k_i_d;
    Kokkos::MDRangePolicy<Kokkos::Rank<3>> xyz_policy({1,1,1},{nz+1,ny+1,nx+1});
    Kokkos::parallel_for("average magnetic field", xyz_policy, KOKKOS_LAMBDA(const int z, const int y, const int x) {
        int zone_x = (x - 1) / nx_zone;
        int zone_y = (y - 1) / ny_zone;
        int zone_z = (z - 1) / nz_zone;
        size_t f0_index = VOXEL(x, y, z, nx, ny, nz);
        float cbx = k_interp(f0_index, interpolator_var::cbx);
        float cby = k_interp(f0_index, interpolator_var::cby);
        float cbz = k_interp(f0_index, interpolator_var::cbz);
        Kokkos::atomic_add(&bfield_avg_d(zone_z, zone_y, zone_x, 0), cbx * incells);
        Kokkos::atomic_add(&bfield_avg_d(zone_z, zone_y, zone_x, 1), cby * incells);
        Kokkos::atomic_add(&bfield_avg_d(zone_z, zone_y, zone_x, 2), cbz * incells);
        });
    Kokkos::deep_copy(bfield_avg_h, bfield_avg_d);

    // Local particle energy spectrum
    typedef Kokkos::View<float****>  View4DType;
    View4DType local_spectrum_d( "local_spectrum", nzones_z, nzones_y, nzones_x, nbins+3 );
    View4DType::HostMirror local_spectrum_h = Kokkos::create_mirror_view( local_spectrum_d );
    auto bfield_sub = Kokkos::subview (local_spectrum_h, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, std::make_pair(0, 3));
    for (int isp = 0; isp < nsp; isp++) {   // loop over species
      species_t *sp = find_species_id(edParams.at(isp)->sp_id, species_list);
      sim_log("computing the distribution function for species "<< sp->name);

      Kokkos::deep_copy(local_spectrum_d, 0.0f);
      const int np = sp->np;
      // Load the particle
      auto& k_particles = sp->k_p_d;
      auto& k_particles_i = sp->k_p_i_d;
      Kokkos::parallel_for("energy_spectrum", Kokkos::RangePolicy < Kokkos::DefaultExecutionSpace > (0, np),
      KOKKOS_LAMBDA (size_t p_index)
      {
        float ux = k_particles(p_index, particle_var::ux);
        float uy = k_particles(p_index, particle_var::uy);
        float uz = k_particles(p_index, particle_var::uz);
        int ii = k_particles_i(p_index);
        float ke = sqrt(1.0 + ux * ux + uy * uy + uz * uz) - 1.0;

        int ebin = floor((log10(ke) - emin_spect_log) / dloge) + 1;
        int zone_z = (ii / nxy2 - 1) / nz_zone;
        int zone_y = ((ii % nxy2) / nx2 - 1) / ny_zone;
        int zone_x = (ii % nx2 - 1) / nx_zone;
        if (ebin < nbins && ebin >= 0) {
          Kokkos::atomic_add(&local_spectrum_d(zone_z, zone_y, zone_x, ebin+3), 1);
        }
      });

      Kokkos::deep_copy(local_spectrum_h, local_spectrum_d);
      Kokkos::deep_copy(bfield_sub, bfield_avg_h); // local magnetic field

      double el2 = uptime();

      sprintf(dset_name, "spectrum_%s", sp->name);
      hid_t dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_FLOAT, filespace,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, local_spectrum_h.data());
      H5Dclose(dset_id);

      el2 = uptime() - el2;
      sim_log("TimeHDF5Write: "<< el2 << " s");

    } // end species loop

    // Finished HDF5 IO
    double el3 = uptime();
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
    el3 = uptime() - el3;
    sim_log("TimeHDF5Close: "<< el3 << " s");
  } // if (should_dump(spectrum))
} // end energy diagnostics
