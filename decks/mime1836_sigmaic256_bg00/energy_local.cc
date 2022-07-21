{
  // energy band diagnostics
  static int energy_diagnostics_first_time = 1;
  static float *dist;         // array to hold the distribution function
  static std::vector<edata *> edParams; // vector to hold species parameters, such as vth
  static std::size_t total_cells;

  // local energy spectrum diagnostics
  static float dloge, emin_spect_log;
  static float *local_spectrum;           // local particle energy spectrum
#ifdef CHECK_PARTICLE_CROSSING_EGTB
  static float *local_spectrum_egtb;      // local particle energy spectrum for particles crossed E>B regions
#endif
  static float *bx_avg, *by_avg, *bz_avg; // local averaged magnetic field
  static int nzones_x;    // number of zones along the x-direction
  static int nzones_y;    // number of zones along the y-direction
  static int nzones_z;    // number of zones along the z-direction
  static int nzones_tot;  // total number of zones in each MPI rank
  static int nbins;       // number of energy bins
  static float incells;   // inverse of total cells in each zone

  const static int nx = grid->nx;
  const static int ny = grid->ny;
  const static int nz = grid->nz;
  const static int nx2 = nx + 2;
  const static int ny2 = ny + 2;
  const static int nz2 = nz + 2;
  const interpolator_t * ALIGNED(16) fi;

  if (should_dump(spectrum)) {
    if (energy_diagnostics_first_time) {
      sim_log("initializing the energy diagnostics");

      double t1 = uptime();

#ifdef ENERGY_BAND_DIAGNOSTICS
      total_cells = nx2 * ny2 * nz2;
      ALLOCATE(dist, global->nbands * total_cells, float);
#endif
      nzones_x = (grid->nx + global->nx_zone - 1) / global->nx_zone;
      nzones_y = (grid->ny + global->ny_zone - 1) / global->ny_zone;
      nzones_z = (grid->nz + global->nz_zone - 1) / global->nz_zone;
      nzones_tot = nzones_x * nzones_y * nzones_z;
      nbins = global->nbins;
      incells = 1.0 / static_cast<float>(global->nx_zone * global->ny_zone * global->nz_zone);

      ALLOCATE(local_spectrum, (nbins+3) * nzones_tot, float);
#ifdef CHECK_PARTICLE_CROSSING_EGTB
      ALLOCATE(local_spectrum_egtb, (nbins+3) * nzones_tot, float);
#endif
      ALLOCATE(bx_avg, nzones_tot, float);
      ALLOCATE(by_avg, nzones_tot, float);
      ALLOCATE(bz_avg, nzones_tot, float);

      emin_spect_log = log10(global->emin_spect);
      dloge = (log10(global->emax_spect) - emin_spect_log) / (nbins - 1);

#ifdef TURBULENCE_MIXING_DIAGNOSTICS
      edParams.push_back(&global->edeTop);
      edParams.push_back(&global->edeBot);
      edParams.push_back(&global->edHTop);
      edParams.push_back(&global->edHBot);
#else
      edParams.push_back(&global->ede);
      edParams.push_back(&global->edH);
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
      edParams.push_back(&global->ede_wo_epara);
      edParams.push_back(&global->edH_wo_epara);
      edParams.push_back(&global->ede_wo_eparay);
      edParams.push_back(&global->edH_wo_eparay);
      edParams.push_back(&global->ede_wo_egtb);
      edParams.push_back(&global->edH_wo_egtb);
#ifdef SPECIES_ADVANCED_WITHOUT_EPERP
      edParams.push_back(&global->ede_wo_eperp);
      edParams.push_back(&global->edH_wo_eperp);
#endif
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
      edParams.push_back(&global->ede_egtb);
      edParams.push_back(&global->edH_egtb);
#endif
#endif

      energy_diagnostics_first_time = 0;

      t1 = uptime() - t1;
      sim_log("Time in initializing energy diagnostics: "<< t1 << " s");
    } // first time called

    // averaged magnetic field in each zone
    CLEAR(bx_avg, nzones_tot);
    CLEAR(by_avg, nzones_tot);
    CLEAR(bz_avg, nzones_tot);
    int zone_x, zone_y, zone_z, zone;
    for (int k = 1; k <= nz; ++k) {
      zone_z = (k - 1) / global->nz_zone;
      for (int j = 1; j <= ny; ++j) {
        zone_y = (j - 1) / global->ny_zone;
        fi = &interpolator(VOXEL(1, j, k, nx, ny, nz));
        for (int i = 1; i <= nx; ++i) {
          zone_x = (i - 1) / global->nx_zone;
          zone = zone_x + nzones_x * (zone_y + zone_z * nzones_y);
          bx_avg[zone] += fi->cbx * incells;
          by_avg[zone] += fi->cby * incells;
          bz_avg[zone] += fi->cbz * incells;
          fi++;
        }
      }
    }

    // Assuming electron and ion are the only two species in species_list
    // Modify if you have more species in the list
    /* h5part_int64_t ierr; */
    /* H5PartFile * h5pf; */
    /* h5part_float32_t *Pf; */
    int npoints_local;

    char fname[256];
    char spectrum_scratch[128];
    char subspectrum_scratch[128];
    char dset_name[64];
    /* sprintf(spectrum_scratch, DUMP_DIR_FORMAT, "spectrum"); */
    sprintf(spectrum_scratch, "spectrum");
    dump_mkdir(spectrum_scratch);
    sprintf(subspectrum_scratch, "%s/T.%d/", spectrum_scratch, step());
    dump_mkdir(subspectrum_scratch);

    int nsp = global->nsp;
    int nsp_pic = nsp;
#ifndef TURBULENCE_MIXING_DIAGNOSTICS
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
    nsp += global->nsp_efield;
#endif
#ifdef CHECK_PARTICLE_CROSSING_EGTB
    nsp += 2;
#endif
#endif

    // Prepare for HDF5 IO
    sprintf(fname, "%s/spectrum_%d.h5part", subspectrum_scratch, step());

    double el1 = uptime();

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(fname , H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    npoints_local = (nbins + 3) * nzones_tot;

    long long total_points, offset;
    long long numpoints = npoints_local;
    MPI_Allreduce(&numpoints, &total_points, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Scan(&numpoints, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    offset -= numpoints;

    hid_t filespace = H5Screate_simple(1, (hsize_t *) &total_points, NULL);
    hid_t memspace =  H5Screate_simple(1, (hsize_t *) &numpoints, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);

    //Comment out for test only
    /* H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); */
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL, (hsize_t *) &numpoints, NULL);
    el1 = uptime() - el1;
    sim_log("TimeHDF5Open: "<< el1 << " s");  //Easy to handle results for scripts

    for (int isp = 0; isp < nsp; isp++) {   // loop over species
      species_t *sp;

      double t2 = uptime();

      if (isp < nsp_pic) {
        sp = find_species_id(edParams.at(isp)->sp_id, species_list);
      } else if (isp < nsp_pic + global->nsp_efield) {
#ifdef SPECIES_ADVANCED_WITH_PART_EFIELD
        sp = find_species_id(edParams.at(isp)->sp_id, global->species_list_efield);
#endif
      } else {
#ifdef CHECK_PARTICLE_CROSSING_EGTB
        sp = find_species_id(edParams.at(isp)->sp_id, global->species_list_egtb);
#endif
      }
      sim_log("computing the distribution function for species "<< sp->name);
      double vth = edParams.at(isp)->vth;
#ifdef ENERGY_BAND_DIAGNOSTICS
      double dke_band = global->emax_band * (vth * vth / 2.0) / global->nbands;
#endif

      CLEAR(local_spectrum, (nbins + 3) * nzones_tot);
#ifdef CHECK_PARTICLE_CROSSING_EGTB
      CLEAR(local_spectrum_egtb, (nbins + 3) * nzones_tot);
#endif
      for (int k = 0; k < nzones_tot; k++) {
        local_spectrum[k*(nbins+3)  ] = bx_avg[k];
        local_spectrum[k*(nbins+3)+1] = by_avg[k];
        local_spectrum[k*(nbins+3)+2] = bz_avg[k];
#ifdef CHECK_PARTICLE_CROSSING_EGTB
        local_spectrum_egtb[k*(nbins+3)  ] = bx_avg[k];
        local_spectrum_egtb[k*(nbins+3)+1] = by_avg[k];
        local_spectrum_egtb[k*(nbins+3)+2] = bz_avg[k];
#endif
      }

#ifdef ENERGY_BAND_DIAGNOSTICS
      CLEAR(dist, total_cells * global->nbands);
#endif
      particle_t *p = sp->p;   // header of the particle array
      int zone_x, zone_y, zone_z, zone;
      int eband;               // energy band
      int ebin;                // energy bin for particle spectrum diagnostics
      for (std::size_t i = 0; i < sp->np; ++i ) { // loop over particles
        double ke = sqrt(1.0 + p->ux*p->ux + p->uy*p->uy + p->uz*p->uz) - 1.0;
#ifdef ENERGY_BAND_DIAGNOSTICS
        eband = int(ke / dke_band);
        // increment the corresponding bin for cell p->i
        if (eband < global->nbands && eband >= 0)
          dist[eband * total_cells + p->i]++;
#endif

        ebin = floor((log10(ke) - emin_spect_log) / dloge) + 1;
        zone_z = (p->i / (nx2 * ny2) - 1) / global->nz_zone;
        zone_y = ((p->i % (nx2 * ny2)) / nx2 - 1) / global->ny_zone;
        zone_x = (p->i % nx2 - 1) / global->nx_zone;
        zone = zone_x + nzones_x * (zone_y + zone_z * nzones_y);
        if (ebin < nbins && ebin >= 0) {
          local_spectrum[zone * (nbins + 3) + ebin + 3]++;
#ifdef CHECK_PARTICLE_CROSSING_EGTB
          if (isp >= (nsp_pic + global->nsp_efield) && p->w > 10) {
            local_spectrum_egtb[zone * (nbins + 3) + ebin + 3]++;
          }
#endif
        }
        p++; // next particle
      }

      t2 = uptime() - t2;
      sim_log("Time in accumulating energy spectrum: "<< t2 << " s");

#ifdef ENERGY_BAND_DIAGNOSTICS

      t3 = uptime():

      // normalize the distribution function
      int ix, iy, iz, ixn, iyn, izn, gcell;

      double np;
      std::size_t icell;
      for (iz = 0; iz < nz2; iz++) {
        for (iy = 0; iy < ny2; iy++) {
          for (ix = 0; ix < nx2; ix++) {
            np = 0;  // particle counter

            icell = LOCAL_CELL_ID(ix, iy, iz);

            for (int k = 0; k < global->nbands; k++) {
              np += dist[k * total_cells + icell];  // count the particles in the cell
            }
            if (np > 0) {
              for ( int k = 0; k < global->nbands; k++) {
                dist[k * total_cells + icell] /= np;  // normalize the distribution
              }
            }

            // is this a ghost cell ?
            gcell = (ix==0) || (ix==nx2-1) || (iy==0) || (iy==ny2-1) || (iz==0) || (iz==nz2-1);

            // if yes, assign the value from a neighboring cell
            if (gcell) {
              // find the neighboring cell
              ixn = ix;
              iyn = iy;
              izn = iz;

              if (ix == 0) ixn++;
              if (ix == nx2-1) ixn--;
              if (iy == 0) iyn++;
              if (iy == ny2-1) iyn--;
              if (iz == 0) izn++;
              if (iz == nz2-1) izn--;

              std::size_t nid = LOCAL_CELL_ID (ixn,iyn,izn);

              for (int k = 0; k < global->nbands; k++) {
                dist[k * total_cells + icell] = dist[k * total_cells + nid];
              }
            } // if (gcell)
          } // ix
        } // iy
      } // iz

      t3 = uptime() - t3;
      sim_log("Time in accumulating hydro for energy bands: "<< t3 << " s");
#endif

#ifdef ENERGY_BAND_DIAGNOSTICS
      // dump all the energy bands
      double t4 = uptime();

      char fname[256];
      sim_log(" writing the distribution function to file ");
      sprintf(fname, HYDRO_FILE_FORMAT, NUMFOLD, step(), edParams.at(isp)->fname, step(), (int)rank());
      sim_log("append data to "<<fname);

      FileIO fedump;
      FileIOStatus status = fedump.open(fname, io_append );
      if (status == fail) ERROR(("Could not open file."));

      fedump.write(dist, total_cells * global->nbands);
      fedump.close();

      t4 = uptime() - t4;
      sim_log("Time in writing hydro for energy bands: "<< t4 << " s");
#endif

      double el2 = uptime();

      sprintf(dset_name, "spectrum_%s", sp->name);
      hid_t dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_FLOAT, filespace,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, local_spectrum);
      H5Dclose(dset_id);

      el2 = uptime() - el2;
      sim_log("TimeHDF5Write: "<< el2 << " s");

#ifdef CHECK_PARTICLE_CROSSING_EGTB
      if (isp >= (nsp_pic + global->nsp_efield)) {
        double el2 = uptime();

        sprintf(dset_name, "spectrum_%s_egtb", sp->name);
        hid_t dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_FLOAT, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, local_spectrum_egtb);
        H5Dclose(dset_id);

        el2 = uptime() - el2;
        sim_log("TimeHDF5Write: "<< el2 << " s");
      }
#endif
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
