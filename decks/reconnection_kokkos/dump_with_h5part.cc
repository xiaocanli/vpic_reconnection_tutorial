/*
 * We should only use H5Part for a small part of the particles because we have
 * to center_p the particles. This requires to make a copy of the part of particle
 * data to be dumped and keep the original data unchanged.
 */
{
    char fname[256], gname[256];
    char particle_scratch[128];
    char subparticle_scratch[128];

    h5part_int64_t ierr;
    int np_local;
    species_t *sp;

    H5PartFile * h5pf;

    h5part_float32_t *Pf ;
    h5part_int32_t *Pi ;

    // get the total number of particles. in this example, output only electrons
    //sp = find_species ("electron");
    sp = species_list;
    /* sprintf(particle_scratch, DUMP_DIR_FORMAT, "particle"); */
    sprintf(particle_scratch, "particle");
    dump_mkdir(particle_scratch);
    sprintf(subparticle_scratch, "%s/T.%d/", particle_scratch, step());
    dump_mkdir(subparticle_scratch);

    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_cb_read", "automatic");
    MPI_Info_set(info, "romio_cb_write", "automatic");

    while (sp) {
        // number of particles on this rank to be dumped
        np_local = (sp->np + global->stride_particle_dump - 1) / global->stride_particle_dump;

        // make a copy of the part of particle data to be dumped
        double ec1 = uptime();

        int sp_np = sp->np;
        int sp_max_np = sp->max_np;
        particle_t * ALIGNED(128) p_buf = NULL;
        if( !p_buf ) MALLOC_ALIGNED( p_buf, np_local, 128 );
        particle_t * sp_p = sp->p;
        sp->p = p_buf;
        sp->np = np_local;
        sp->max_np = np_local;

        for (long long iptl = 0, i = 0; iptl < sp_np; iptl += global->stride_particle_dump, ++i) {
            COPY( &sp->p[i], &sp_p[iptl], 1 );
        }
        center_p(sp, interpolator_array);
        ec1 = uptime() - ec1;
        sim_log("Time in copying particle data: "<< ec1 << " s");

        Pf = (h5part_float32_t *) sp->p;
        Pi = (h5part_int32_t *) sp->p;

        // open H5part file in "particle/T.<step>/" subdirectory
        // filename: eparticle.h5p
        sprintf (fname, "%s/%s_%d.h5part", subparticle_scratch, sp->name, step());

        double el1 = uptime();
        // h5pf = H5PartOpenFileParallel (fname, H5PART_WRITE | H5PART_FS_LUSTRE, MPI_COMM_WORLD);
        h5pf = H5PartOpenFileParallel (fname, H5PART_WRITE, MPI_COMM_WORLD);
        // h5pf = H5PartOpenFileParallel (fname, H5PART_WRITE | H5PART_VFD_MPIIO_IND, MPI_COMM_WORLD);

        ierr = H5PartSetStep (h5pf, step());
        ierr = H5PartSetNumParticlesStrided (h5pf, np_local, 8);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in setting the # of particles; Error code: " << ierr);

        el1 = uptime() - el1;
        sim_log("Time in opening H5part file for particle dump: "<< el1 << " s");

        double el2 = uptime();

        ierr = H5PartWriteDataFloat32 (h5pf, "dX", Pf);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dX; Error code: " << ierr);
        sim_log("Finished writing dX");

        ierr = H5PartWriteDataFloat32 (h5pf, "dY", Pf+1);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dY; Error code: " << ierr);
        sim_log("Finished writing dY");

        ierr = H5PartWriteDataFloat32 (h5pf, "dZ", Pf+2);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dZ; Error code: " << ierr);
        sim_log("Finished writing dZ");

        ierr = H5PartWriteDataInt32   (h5pf, "i",  Pi+3);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing i; Error code: " << ierr);
        sim_log("Finished writing i");

        ierr = H5PartWriteDataFloat32 (h5pf, "Ux", Pf+4);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Ux; Error code: " << ierr);
        sim_log("Finished writing Ux");

        ierr = H5PartWriteDataFloat32 (h5pf, "Uy", Pf+5);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Uy; Error code: " << ierr);
        sim_log("Finished writing Uy");

        ierr = H5PartWriteDataFloat32 (h5pf, "Uz", Pf+6);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing Uz; Error code: " << ierr);
        sim_log("Finished writing Uz");

        ierr = H5PartWriteDataFloat32 (h5pf, "q", Pf+7);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing q; Error code: " << ierr);
        sim_log("Finished writing q");

        el2 = uptime() - el2;
        sim_log("Time in writing H5part file for particle dump: "<< el2 << " s");

        double el3 = uptime();
        ierr = H5PartCloseFile (h5pf);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in closing " << fname << ". Error code: " << ierr);

        el3 = uptime() - el3;
        sim_log("Time in closing H5part file for particle dump: "<< el3 << " s");

        sp->p      = sp_p;
        sp->np     = sp_np;
        sp->max_np = sp_max_np;
        FREE_ALIGNED(p_buf);

        // Write metadata if step() == 0
        char meta_fname[256];

        sprintf (meta_fname, "%s/grid_metadata_%s_%d.h5part", subparticle_scratch,
                sp->name, step());

#if 0
        double el4 = uptime();
        H5PartFile * h5pmf;

        h5pmf = H5PartOpenFileParallel (meta_fname, H5PART_WRITE, MPI_COMM_WORLD);
        // h5pmf = H5PartOpenFileParallel (meta_fname, H5PART_WRITE | H5PART_VFD_MPIIO_IND, MPI_COMM_WORLD);

        ierr = H5PartSetStep (h5pmf, step());
        ierr = H5PartSetNumParticles (h5pmf, 1);
        if (ierr != H5PART_SUCCESS) {
            sim_log ("Error occured in setting the # of grid metadata entries; Error code: " << ierr);
        }

        ierr = H5PartWriteDataInt32(h5pmf, "np_local", (h5part_int32_t *) &np_local);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing np_local; Error code: " << ierr);

        ierr = H5PartWriteDataInt32(h5pmf, "nx", (h5part_int32_t *) &grid->nx);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing nx; Error code: " << ierr);

        ierr = H5PartWriteDataInt32(h5pmf, "ny", (h5part_int32_t *) &grid->ny);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing ny; Error code: " << ierr);

        ierr = H5PartWriteDataInt32(h5pmf, "nz", (h5part_int32_t* ) &grid->nz);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing nz; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "x0", (h5part_float32_t* ) &grid->x0);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing x0; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "y0", (h5part_float32_t* ) &grid->y0);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing y0; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "z0", (h5part_float32_t* ) &grid->z0);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing z0; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "dx", (h5part_float32_t* ) &grid->dx);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dx; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "dy", (h5part_float32_t* ) &grid->dy);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dy; Error code: " << ierr);

        ierr = H5PartWriteDataFloat32(h5pmf, "dz", (h5part_float32_t* ) &grid->dz);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in writing dz; Error code: " << ierr);

        ierr = H5PartCloseFile (h5pmf);
        if (ierr != H5PART_SUCCESS) sim_log ("Error occured in closing " << meta_fname << ". Error code: " << ierr);
        el4 = uptime() - el4;
        sim_log("Time in writing metadata: "<< el4 << " s");
#endif
        sprintf(gname, "Step#%d", step());

        double el4 = uptime();

        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
        hid_t file_id = H5Fcreate(meta_fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        hid_t group_id = H5Gcreate(file_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Pclose(plist_id);

        hsize_t mpi_rank, mpi_size, ndata;
        mpi_rank = rank();
        mpi_size = nproc();
        ndata = 1;
        hid_t filespace = H5Screate_simple(1, &mpi_size, NULL);
        hid_t memspace =  H5Screate_simple(1, &ndata, NULL);

        plist_id = H5Pcreate(H5P_DATASET_XFER);
        /* H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); */
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &mpi_rank, NULL, &ndata, NULL);

        el4 = uptime() - el4;
        sim_log("Time in opening HDF5 file for particle metafile: "<< el4 << " s");

        double el5 = uptime();

        hid_t dset_id = H5Dcreate(group_id, "np_local", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &np_local);
        H5Dclose(dset_id);
        /* sim_log("Finished writing np_local"); */

        dset_id = H5Dcreate(group_id, "nx", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &grid->nx);
        H5Dclose(dset_id);
        /* sim_log("Finished writing nx"); */

        dset_id = H5Dcreate(group_id, "ny", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &grid->ny);
        H5Dclose(dset_id);
        /* sim_log("Finished writing ny"); */

        dset_id = H5Dcreate(group_id, "nz", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &grid->nz);
        H5Dclose(dset_id);
        /* sim_log("Finished writing nz"); */

        dset_id = H5Dcreate(group_id, "x0", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->x0);
        H5Dclose(dset_id);
        /* sim_log("Finished writing x0"); */

        dset_id = H5Dcreate(group_id, "y0", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->y0);
        H5Dclose(dset_id);
        /* sim_log("Finished writing y0"); */

        dset_id = H5Dcreate(group_id, "z0", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->z0);
        H5Dclose(dset_id);
        /* sim_log("Finished writing z0"); */

        dset_id = H5Dcreate(group_id, "dx", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->dx);
        H5Dclose(dset_id);
        /* sim_log("Finished writing dx"); */

        dset_id = H5Dcreate(group_id, "dy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->dy);
        H5Dclose(dset_id);
        /* sim_log("Finished writing dy"); */

        dset_id = H5Dcreate(group_id, "dz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, &grid->dz);
        H5Dclose(dset_id);
        /* sim_log("Finished writing dz"); */

        el5 = uptime() - el5;
        sim_log("Time in writing HDF5 file for particle metafile: "<< el5 << " s");

        double el6 = uptime();

        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Pclose(plist_id);
        H5Gclose(group_id);
        H5Fclose(file_id);

        el6 = uptime() - el6;
        sim_log("Time in closing HDF5 file for particle metafile: "<< el6 << " s");

        sp = sp->next;
    }

    MPI_Info_free(&info);
}
