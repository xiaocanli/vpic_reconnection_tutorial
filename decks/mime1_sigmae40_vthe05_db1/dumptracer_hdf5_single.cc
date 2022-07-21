{
  const int nx2 = grid->nx + 2;
  const int nx2ny2 = (grid->ny+2) * nx2;
  const particle_t     * ALIGNED(32) p;
  const interpolator_t * ALIGNED(16) f0=interpolator_array->i;
  const interpolator_t * ALIGNED(16) f;
  const grid_t * g = grid;
  const hydro_t * ALIGNED(32) hy;
  const float r8V = 0.125;
  float ex, ey, ez, bx, by, bz; // EMF at particle position
  float vx, vy, vz; // Single-fluid velocity at particle position
  float vex, vey, vez; // Electron bulk velocity
  float ne, ni; // Number densities at particle position
  float dx0, dy0, dz0;
  float dx, dy, dz;
  float w0, w1, w2, w3, w4, w5, w6, w7;
  int ii, j;
  int ix, iy, iz;
  int hydro_data_shift;

  if (global->emf_at_tracer) {
    hydro_data_shift = 14;
  } else {
    hydro_data_shift = 8;
  }

  // prepare the tracer particle data
  static hid_t file_id = -1;
  static std::list<std::list<std::vector<float>>> tracers;
  static std::list<std::vector<int>> ntracer_particles;
  species_t *sp = global->tracers_list;
  std::size_t nvar = 8 + global->num_tracer_fields_add;
  if (file_id < 0) {
    tracers.clear();
    ntracer_particles.clear();
    while (sp) {
      tracers.push_back(std::list<std::vector<float>>());
      for (std::size_t i = 0; i != nvar; ++i) {
        tracers.back().push_back(std::vector<float>());
      }
      ntracer_particles.push_back(std::vector<int>());
      sp = sp->next;
    }
    file_id = 1;
  }
  static std::vector<std::string> var_names = {"dX", "dY", "dZ", "i", "Ux", "Uy", "Uz", "q"};
  if ( global->emf_at_tracer ) {
    var_names.push_back("Ex");
    var_names.push_back("Ey");
    var_names.push_back("Ez");
    var_names.push_back("Bx");
    var_names.push_back("By");
    var_names.push_back("Bz");
  }
  if ( global->hydro_at_tracer ) {
    var_names.push_back("Vx");
    var_names.push_back("Vy");
    var_names.push_back("Vz");
    var_names.push_back("ne");
    var_names.push_back("ni");
    if ( global->ve_at_tracer ) {
      var_names.push_back("Vex");
      var_names.push_back("Vey");
      var_names.push_back("Vez");
    }
  }

  sp = global->tracers_list;
  auto tracer_species = tracers.begin();
  auto nparticle = ntracer_particles.begin();
  while (sp) {
    int np_local = sp->np; // number of particles on this rank
    float *pout = new float[np_local * nvar];

    // Gather particle information
    for (p = sp->p, j = 0; j < np_local; j++, p++ ) {
      dx0 = p->dx;
      dy0 = p->dy;
      dz0 = p->dz;
      ii  = p->i;
      f   = f0 + ii;
      CALC_TRACER_USER_DEFINED_DATA;
      pout[j]  = tracer_x;
      pout[j + np_local * 1]  = tracer_y;
      pout[j + np_local * 2]  = tracer_z;
      pout[j + np_local * 3]  = *reinterpret_cast<float*>(&ii);
      pout[j + np_local * 4]  = p->ux;
      pout[j + np_local * 5]  = p->uy;
      pout[j + np_local * 6]  = p->uz;
      pout[j + np_local * 7]  = p->w;
      if ( global->emf_at_tracer ) {
        pout[j + np_local * 8]  = ex;
        pout[j + np_local * 9]  = ey;
        pout[j + np_local * 10]  = ez;
        pout[j + np_local * 11]  = bx;
        pout[j + np_local * 12]  = by;
        pout[j + np_local * 13]  = bz;
        if ( global->hydro_at_tracer ) {
          pout[j + np_local * 14] = vx;
          pout[j + np_local * 15] = vy;
          pout[j + np_local * 16] = vz;
          pout[j + np_local * 17] = ne;
          pout[j + np_local * 18] = ni;
          if ( global->ve_at_tracer ) {
            pout[j + np_local * 19] = vex;
            pout[j + np_local * 20] = vey;
            pout[j + np_local * 21] = vez;
          }
        }
      } else {
        if ( global->hydro_at_tracer ) {
          pout[j + np_local * 8] = vx;
          pout[j + np_local * 9] = vy;
          pout[j + np_local * 10] = vz;
          pout[j + np_local * 11] = ne;
          pout[j + np_local * 12] = ni;
          if ( global->ve_at_tracer ) {
            pout[j + np_local * 13] = vex;
            pout[j + np_local * 14] = vey;
            pout[j + np_local * 15] = vez;
          }
        }
      }
    }

    (*nparticle).push_back(np_local);
    int ivar = 0;
    for (auto tracer = (*tracer_species).begin(); tracer != (*tracer_species).end(); ++tracer) {
      (*tracer).insert((*tracer).end(), pout+ivar*np_local, pout+(ivar+1)*np_local);
      ++ivar;
    }
    delete [] pout;

    sp = sp->next;
    ++nparticle;
    ++tracer_species;
  }

  // whether simulation reaches quota time
  bool reach_quota_time = step()>0 && global->quota_check_interval &&
    (step()%global->quota_check_interval) == 0 && uptime() > global->quota_sec;
  bool to_dump_tracer = (step() + global->tracer_interval) % global->tracer_file_interval == 0;
  if (step() == num_step || to_dump_tracer || reach_quota_time) {
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_cb_read", "automatic");
    MPI_Info_set(info, "romio_cb_write", "automatic");

    double el1 = uptime();

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    MPI_Info_free(&info);

    char fname[256], gname[256];
    char tracer_dir[64];
    //sprintf(tracer_dir, "tracer/tracer%d/T.%d", global->particle_tracing, step());
    sprintf(tracer_dir, "tracer/tracer1/T.%d",
        (step() / global->tracer_file_interval) * global->tracer_file_interval);
    dump_mkdir(tracer_dir);
    sprintf(fname, "%s/tracers.h5p", tracer_dir);
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    sprintf(gname, "Step#%d", (step() / global->tracer_file_interval) * global->tracer_file_interval);
    hid_t group_id = H5Gcreate(file_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t meta_group_id = H5Gcreate(group_id, "grid_metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    el1 = uptime() - el1;
    sim_log("Time in opening HDF5 file for tracers: "<< el1 << " s");

    double el2 = uptime();

    auto tracer_species = tracers.begin();
    auto nparticle = ntracer_particles.begin();
    sp = global->tracers_list;
    while (sp) {
      hid_t sub_group_id = H5Gcreate(group_id, sp->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      long long total_particles, offset;
      long long numparticles = (*tracer_species).front().size();
      MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      offset -= numparticles;

      plist_id = H5Pcreate(H5P_DATASET_XFER);
      hid_t filespace, memspace;

      if (total_particles > 0) {
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        /* H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */

        filespace = H5Screate_simple(1, (hsize_t *) &total_particles, NULL);
        memspace =  H5Screate_simple(1, (hsize_t *) &numparticles, NULL);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL,
            (hsize_t *) &numparticles, NULL);

        int ivar = 0;
        for (auto tracer = (*tracer_species).begin(); tracer != (*tracer_species).end(); ++tracer) {
          float *Pf = (float *) ((*tracer).data());
          int *Pi = (int *) ((*tracer).data());
          std::string var_name = var_names[ivar];

          hid_t dset_id;
          if (var_name == "i" || var_name == "q") {
            dset_id = H5Dcreate(sub_group_id, var_name.c_str(), H5T_NATIVE_INT,
                filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            int ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, Pi);
          } else {
            dset_id = H5Dcreate(sub_group_id, var_name.c_str(), H5T_NATIVE_FLOAT,
                filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf);
          }
          H5Dclose(dset_id);

          ++ivar;
        }

        H5Sclose(memspace);
        H5Sclose(filespace);
      }

      H5Gclose(sub_group_id);

      // Meta data
      hsize_t dcount[2], doffset[2], dset_dims[2];
      dcount[0] = (*nparticle).size();
      dcount[1] = 1;
      doffset[0] = 0;
      doffset[1] = rank();
      dset_dims[0] = dcount[0];
      dset_dims[1] = nproc();

      filespace = H5Screate_simple(2, dset_dims, NULL);
      memspace =  H5Screate_simple(2, dcount, NULL);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, doffset, NULL, dcount, NULL);

      char meta_name[64];
      sprintf(meta_name, "np_local_%s", sp->name);
      hid_t dset_id = H5Dcreate(meta_group_id, meta_name, H5T_NATIVE_INT, filespace,
              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      int ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
          plist_id, (*nparticle).data());
      H5Dclose(dset_id);

      H5Sclose(memspace);
      H5Sclose(filespace);

      H5Pclose(plist_id);
      sp = sp->next;
      ++nparticle;
      ++tracer_species;
    }

    el2 = uptime() - el2;
    sim_log("Time in writing HDF5 file for tracers: "<< el2 << " s");

    double el3 = uptime();

    H5Gclose(meta_group_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    el3 = uptime() - el3;
    sim_log("Time in closing HDF5 file for tracers: "<< el3 << " s");

    file_id = -1;
  }

}
