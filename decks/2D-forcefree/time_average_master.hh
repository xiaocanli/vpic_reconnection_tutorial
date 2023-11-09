/*--------------------------------------------------------------
 * time averaging header file for VPIC master version
 *--------------------------------------------------------------*/
#ifndef __TIME_AVERAGE_MASTER_HH__
#define __TIME_AVERAGE_MASTER_HH__

// hydro at current step
#ifndef HYDRO_STEP
#define HYDRO_STEP(x,y,z) hi[LOCAL_CELL_ID(x,y,z)]
#endif

// fields at current step
#ifndef FIELD_STEP
#define FIELD_STEP(x,y,z) fi[LOCAL_CELL_ID(x,y,z)]
#endif

// averaged electron hydro at current step
#ifndef HYDRO_AVG_E
#define HYDRO_AVG_E(x,y,z) hydro_avg_e[LOCAL_CELL_ID(x,y,z)]
#endif

#ifndef HYDRO_AVG_ETOP
#define HYDRO_AVG_ETOP(x,y,z) hydro_avg_etop[LOCAL_CELL_ID(x,y,z)]
#endif

#ifndef HYDRO_AVG_EBOT
#define HYDRO_AVG_EBOT(x,y,z) hydro_avg_ebot[LOCAL_CELL_ID(x,y,z)]
#endif

// averaged ion hydro at current step
#ifndef HYDRO_AVG_H
#define HYDRO_AVG_H(x,y,z) hydro_avg_h[LOCAL_CELL_ID(x,y,z)]
#endif

#ifndef HYDRO_AVG_HTOP
#define HYDRO_AVG_HTOP(x,y,z) hydro_avg_htop[LOCAL_CELL_ID(x,y,z)]
#endif

#ifndef HYDRO_AVG_HBOT
#define HYDRO_AVG_HBOT(x,y,z) hydro_avg_hbot[LOCAL_CELL_ID(x,y,z)]
#endif

// averaged fields at current step
#ifndef EMF_AVG
#define EMF_AVG(x,y,z) emf_avg[LOCAL_CELL_ID(x,y,z)]
#endif

// loop over yz face
#ifndef LOOP_YZ_FACE
#define LOOP_YZ_FACE             \
  for (z = 0; z < nz + 2; z++)   \
    for (y = 0; y < ny + 2; y++)
#endif

// update hydro along x
#ifndef HYDRO_ALONG_X
#define HYDRO_ALONG_X            \
  for (x = 0; x < nx + 2; x++) { \
    havg->jx += h0->jx;          \
    havg->jy += h0->jy;          \
    havg->jz += h0->jz;          \
    havg->rho += h0->rho;        \
    havg->px += h0->px;          \
    havg->py += h0->py;          \
    havg->pz += h0->pz;          \
    havg->ke += h0->ke;          \
    havg->txx += h0->txx;        \
    havg->tyy += h0->tyy;        \
    havg->tzz += h0->tzz;        \
    havg->tyz += h0->tyz;        \
    havg->tzx += h0->tzx;        \
    havg->txy += h0->txy;        \
    havg++;                      \
    h0++;                        \
  }
#endif

// update fields along x
#ifndef EMF_ALONG_X
#define EMF_ALONG_X              \
  for (x = 0; x < nx + 2; x++) { \
    favg->ex += f0->ex;          \
    favg->ey += f0->ey;          \
    favg->ez += f0->ez;          \
    favg->cbx += f0->cbx;        \
    favg->cby += f0->cby;        \
    favg->cbz += f0->cbz;        \
    favg++;                      \
    f0++;                        \
  }
#endif

// calculate the averaged hydro
#ifndef AVG_HYDRO
#define AVG_HYDRO                \
  for (x = 0; x < nx + 2; x++) { \
    havg->jx *= isteps;          \
    havg->jy *= isteps;          \
    havg->jz *= isteps;          \
    havg->rho *= isteps;         \
    havg->px *= isteps;          \
    havg->py *= isteps;          \
    havg->pz *= isteps;          \
    havg->ke *= isteps;          \
    havg->txx *= isteps;         \
    havg->tyy *= isteps;         \
    havg->tzz *= isteps;         \
    havg->tyz *= isteps;         \
    havg->tzx *= isteps;         \
    havg->txy *= isteps;         \
    havg++;                      \
  }
#endif

// calculate the averaged fields
#ifndef AVG_EMF
#define AVG_EMF                  \
  for (x = 0; x < nx + 2; x++) { \
    favg->ex *= isteps;          \
    favg->ey *= isteps;          \
    favg->ez *= isteps;          \
    favg->cbx *= isteps;         \
    favg->cby *= isteps;         \
    favg->cbz *= isteps;         \
    favg++;                      \
  }
#endif

// dump hydro in binary format (N-to-N)
#ifndef DUMP_HYDRO_AVG
#define DUMP_HYDRO_AVG(filename, species, sname, ha)                              \
  sprintf(fdir, "%s/T.%d", global->hydro_dir, step_dump);                         \
  dump_mkdir(fdir);                                                               \
  sprintf(fname, "%s/%s.%ld.%d", fdir, filename, step_dump, rank());              \
  status = fileIO.open(fname, io_write);                                          \
  if(status == fail) ERROR(("Failed opening file: %s", fname));                   \
  sp = find_species_name(species, species_list);                                  \
  WRITE_HEADER_V0(2, sp->id, sp->q/sp->m, step(), fileIO);                        \
  dim[0] = nx + 2;                                                                \
  dim[1] = ny + 2;                                                                \
  dim[2] = nz + 2;                                                                \
  WRITE_ARRAY_HEADER(hydro_array->h, 3, dim, fileIO);                             \
  /* Create a variable list of hydro values to output */                          \
  numvars = std::min(global->h##sname##dParams.output_vars.bitsum(),              \
                            total_hydro_variables);                               \
  varlist = new size_t[numvars];                                                  \
  for(size_t i(0), c(0); i<total_hydro_variables; i++)                            \
    if( global->h##sname##dParams.output_vars.bitset(i) ) varlist[c++] = i;       \
                                                                                  \
  for(size_t v(0); v<numvars; v++)                                                \
  for(size_t k(0); k<nz+2; k++)                                                   \
  for(size_t j(0); j<ny+2; j++)                                                   \
  for(size_t i(0); i<nx+2; i++) {                                                 \
    const uint32_t * href = reinterpret_cast<uint32_t *>(&ha(i,j,k));             \
    fileIO.write(&href[varlist[v]], 1);                                           \
  }                                                                               \
                                                                                  \
  if( fileIO.close() ) ERROR(( "File close failed on hydro dump!!!" ));           \
  delete[] varlist;
#endif

// dump fields in binary format (N-to-N)
#ifndef DUMP_EMF_AVG
#define DUMP_EMF_AVG(filename, fa)                                                \
  sprintf(fdir, "%s/T.%d", global->fields_dir, step_dump);                        \
  dump_mkdir(fdir);                                                               \
  sprintf(fname, "%s/%s.%ld.%d", fdir, filename, step_dump, rank());              \
  status = fileIO.open(fname, io_write);                                          \
  if(status == fail) ERROR(("Failed opening file: %s", fname));                   \
  WRITE_HEADER_V0(1, -1, 0, step(), fileIO);                                      \
  dim[0] = nx + 2;                                                                \
  dim[1] = ny + 2;                                                                \
  dim[2] = nz + 2;                                                                \
  WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);                             \
  for(size_t v(0); v<numvars; v++)                                                \
  for(size_t k(0); k<nz+2; k++)                                                   \
  for(size_t j(0); j<ny+2; j++)                                                   \
  for(size_t i(0); i<nx+2; i++) {                                                 \
    const uint32_t * fref = reinterpret_cast<uint32_t *>(&fa(i,j,k));             \
    fileIO.write(&fref[v], 1);                                                    \
  }                                                                               \
                                                                                  \
  if( fileIO.close() ) ERROR(( "File close failed on field dump!!!" ));
#endif

// access to both fields and hydro data
#ifndef FIELD_HYDRO_ARRAY
#define FIELD_HYDRO_ARRAY(DATA_STRUCT, x, y, z) DATA_STRUCT[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]
#endif

// dump fields or hydro data in HDF5 format (N-to-1)
#ifndef DUMP_TO_HDF5
#define DUMP_TO_HDF5(DATA_STRUCT, DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)         \
{                                                                                  \
  dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace,                \
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);                                      \
  temp_buf_index = 0;                                                              \
  for (size_t i(1); i < grid->nx + 1; i++) {                                       \
    for (size_t j(1); j < grid->ny + 1; j++) {                                     \
      for (size_t k(1); k < grid->nz + 1; k++) {                                   \
        temp_buf[temp_buf_index] =                                                 \
            FIELD_HYDRO_ARRAY(DATA_STRUCT, i, j, k).ATTRIBUTE_NAME;                \
        temp_buf_index = temp_buf_index + 1;                                       \
      }                                                                            \
    }                                                                              \
  }                                                                                \
  dataspace_id = H5Dget_space(dset_id);                                            \
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL,                  \
      local_size, NULL);                                                           \
  H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id,                \
      temp_buf);                                                                   \
  H5Sclose(dataspace_id);                                                          \
  H5Dclose(dset_id);                                                               \
}
#endif

// create a HDF5 file in parallel
#ifndef CREATE_HDF5_FILE
#define CREATE_HDF5_FILE(fbase, filename, step_dump)                                \
{                                                                                   \
    char fname[256];                                                                \
    char output_dir[128];                                                           \
                                                                                    \
    sprintf(output_dir, "%s/T.%zu/", fbase, step_dump);                             \
    dump_mkdir(output_dir);                                                         \
                                                                                    \
    sprintf(fname, "%s/%s_%zu.h5", output_dir, filename, step_dump);                \
    plist_id = H5Pcreate(H5P_FILE_ACCESS);                                          \
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);                      \
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);               \
    H5Pclose(plist_id);                                                             \
                                                                                    \
    sprintf(fname, "Timestep_%zu", step_dump);                                      \
    group_id = H5Gcreate(file_id, fname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    \
                                                                                    \
    plist_id = H5Pcreate(H5P_DATASET_XFER);                                         \
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);                               \
                                                                                    \
    global_size[0] = grid->nx * grid->gpx;                                          \
    global_size[1] = grid->ny * grid->gpy;                                          \
    global_size[2] = grid->nz * grid->gpz;                                          \
                                                                                    \
    local_size[0] = grid->nx;                                                       \
    local_size[1] = grid->ny;                                                       \
    local_size[2] = grid->nz;                                                       \
                                                                                    \
    int mpi_rank_x, mpi_rank_y, mpi_rank_z;                                         \
    mpi_rank_x = rank();                                                            \
    mpi_rank_y  = mpi_rank_x / grid->gpx;                                           \
    mpi_rank_x -= mpi_rank_y * grid->gpx;                                           \
    mpi_rank_z  = mpi_rank_y / grid->gpy;                                           \
    mpi_rank_y -= mpi_rank_z * grid->gpy;                                           \
                                                                                    \
    offset[0] = grid->nx * mpi_rank_x;                                              \
    offset[1] = grid->ny * mpi_rank_y;                                              \
    offset[2] = grid->nz * mpi_rank_z;                                              \
                                                                                    \
    filespace = H5Screate_simple(3, global_size, NULL);                             \
    memspace = H5Screate_simple(3, local_size, NULL);                               \
    dataspace_id;                                                                   \
}
#endif

// dump time-averaged fields in HDF5
#ifndef DUMP_TIME_AVERAGE_FIELDS_HDF5
#define DUMP_TIME_AVERAGE_FIELDS_HDF5(fbase, DATA_STRUCT, step_dump)               \
{                                                                                  \
    hid_t plist_id, file_id, group_id, dset_id;                                    \
    hid_t filespace, memspace, dataspace_id;                                       \
    hsize_t temp_buf_index;                                                        \
    hsize_t global_size[3], local_size[3], offset[3];                              \
    float *temp_buf = (float *)malloc(sizeof(float) *                              \
        (grid->nx) * (grid->ny) * (grid->nz));                                     \
    double el2 = uptime();                                                         \
    CREATE_HDF5_FILE(fbase, "fields", step_dump)                                   \
    if (field_dump_flag.ex)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "ex", ex, H5T_NATIVE_FLOAT);                     \
    if (field_dump_flag.ey)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "ey", ey, H5T_NATIVE_FLOAT);                     \
    if (field_dump_flag.ez)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "ez", ez, H5T_NATIVE_FLOAT);                     \
                                                                                   \
    if (field_dump_flag.cbx)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "cbx", cbx, H5T_NATIVE_FLOAT);                   \
    if (field_dump_flag.cby)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "cby", cby, H5T_NATIVE_FLOAT);                   \
    if (field_dump_flag.cbz)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "cbz", cbz, H5T_NATIVE_FLOAT);                   \
                                                                                   \
    el2 = uptime() - el2;                                                          \
    sim_log("TimeHDF5Write: " << el2 << " s");                                     \
                                                                                   \
    double el3 = uptime();                                                         \
                                                                                   \
    free(temp_buf);                                                                \
    H5Sclose(filespace);                                                           \
    H5Sclose(memspace);                                                            \
    H5Pclose(plist_id);                                                            \
    H5Gclose(group_id);                                                            \
    H5Fclose(file_id);                                                             \
                                                                                   \
    el3 = uptime() - el3;                                                          \
    sim_log("TimeHDF5Close: " << el3 << " s");                                     \
}
#endif

// dump time-averaged hydro in HDF5
#ifndef DUMP_TIME_AVERAGE_HYDRO_HDF5
#define DUMP_TIME_AVERAGE_HYDRO_HDF5(speciesname, fbase, DATA_STRUCT, step_dump)   \
{                                                                                  \
    species_t *sp = find_species_name(speciesname, species_list);                  \
    if (!sp)                                                                       \
      ERROR(("Invalid species name: %s", speciesname));                            \
                                                                                   \
    char hname[32];                                                                \
    sprintf(hname, "hydro_%s", speciesname);                                       \
    hid_t plist_id, file_id, group_id, dset_id;                                    \
    hid_t filespace, memspace, dataspace_id;                                       \
    hsize_t temp_buf_index;                                                        \
    hsize_t global_size[3], local_size[3], offset[3];                              \
    float *temp_buf = (float *)malloc(sizeof(float) *                              \
        (grid->nx) * (grid->ny) * (grid->nz));                                     \
    double el2 = uptime();                                                         \
    CREATE_HDF5_FILE(fbase, hname, step_dump)                                      \
                                                                                   \
    if (hydro_dump_flag.jx)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "jx", jx, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.jy)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "jy", jy, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.jz)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "jz", jz, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.rho)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "rho", rho, H5T_NATIVE_FLOAT);                   \
                                                                                   \
    if (hydro_dump_flag.px)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "px", px, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.py)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "py", py, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.pz)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "pz", pz, H5T_NATIVE_FLOAT);                     \
    if (hydro_dump_flag.ke)                                                        \
        DUMP_TO_HDF5(DATA_STRUCT, "ke", ke, H5T_NATIVE_FLOAT);                     \
                                                                                   \
    if (hydro_dump_flag.txx)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "txx", txx, H5T_NATIVE_FLOAT);                   \
    if (hydro_dump_flag.tyy)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "tyy", tyy, H5T_NATIVE_FLOAT);                   \
    if (hydro_dump_flag.tzz)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "tzz", tzz, H5T_NATIVE_FLOAT);                   \
                                                                                   \
    if (hydro_dump_flag.tyz)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "tyz", tyz, H5T_NATIVE_FLOAT);                   \
    if (hydro_dump_flag.tzx)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "tzx", tzx, H5T_NATIVE_FLOAT);                   \
    if (hydro_dump_flag.txy)                                                       \
        DUMP_TO_HDF5(DATA_STRUCT, "txy", txy, H5T_NATIVE_FLOAT);                   \
                                                                                   \
    el2 = uptime() - el2;                                                          \
    sim_log("TimeHDF5Write: " << el2 << " s");                                     \
                                                                                   \
    double el3 = uptime();                                                         \
                                                                                   \
    free(temp_buf);                                                                \
    H5Sclose(filespace);                                                           \
    H5Sclose(memspace);                                                            \
    H5Pclose(plist_id);                                                            \
    H5Gclose(group_id);                                                            \
    H5Fclose(file_id);                                                             \
                                                                                   \
    el3 = uptime() - el3;                                                          \
    sim_log("TimeHDF5Close: " << el3 << " s");                                     \
}
#endif

// read fields or hydro data in HDF5 format
#ifndef READ_FROM_HDF5
#define READ_FROM_HDF5(DATA_STRUCT, DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)       \
{                                                                                  \
  dset_id = H5Dopen(group_id, DSET_NAME, H5P_DEFAULT);                             \
  temp_buf_index = 0;                                                              \
  dataspace_id = H5Dget_space(dset_id);                                            \
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL,                  \
      local_size, NULL);                                                           \
  H5Dread(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id, temp_buf);      \
  H5Sclose(dataspace_id);                                                          \
  H5Dclose(dset_id);                                                               \
  for (size_t i(1); i < grid->nx + 1; i++) {                                       \
    for (size_t j(1); j < grid->ny + 1; j++) {                                     \
      for (size_t k(1); k < grid->nz + 1; k++) {                                   \
        FIELD_HYDRO_ARRAY(DATA_STRUCT, i, j, k).ATTRIBUTE_NAME =                   \
            temp_buf[temp_buf_index];                                              \
        temp_buf_index = temp_buf_index + 1;                                       \
      }                                                                            \
    }                                                                              \
  }                                                                                \
}
#endif

// open a HDF5 file in parallel
#ifndef OPEN_HDF5_FILE
#define OPEN_HDF5_FILE(fbase, filename, time_step)                                  \
{                                                                                   \
    char fname[256];                                                                \
    char reastart_avg_dir[128];                                                     \
                                                                                    \
    sprintf(reastart_avg_dir, "%s/T.%zu/", fbase, time_step);                       \
    dump_mkdir(reastart_avg_dir);                                                   \
                                                                                    \
    sprintf(fname, "%s/%s_%zu.h5", reastart_avg_dir, filename, time_step);          \
    plist_id = H5Pcreate(H5P_FILE_ACCESS);                                          \
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);                      \
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);                             \
    H5Pclose(plist_id);                                                             \
                                                                                    \
    sprintf(fname, "Timestep_%zu", time_step);                                      \
    group_id = H5Gopen(file_id, fname, H5P_DEFAULT);                                \
                                                                                    \
    plist_id = H5Pcreate(H5P_DATASET_XFER);                                         \
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);                               \
                                                                                    \
    global_size[0] = grid->nx * grid->gpx;                                          \
    global_size[1] = grid->ny * grid->gpy;                                          \
    global_size[2] = grid->nz * grid->gpz;                                          \
                                                                                    \
    local_size[0] = grid->nx;                                                       \
    local_size[1] = grid->ny;                                                       \
    local_size[2] = grid->nz;                                                       \
                                                                                    \
    int mpi_rank_x, mpi_rank_y, mpi_rank_z;                                         \
    mpi_rank_x = rank();                                                            \
    mpi_rank_y  = mpi_rank_x / grid->gpx;                                           \
    mpi_rank_x -= mpi_rank_y * grid->gpx;                                           \
    mpi_rank_z  = mpi_rank_y / grid->gpy;                                           \
    mpi_rank_y -= mpi_rank_z * grid->gpy;                                           \
                                                                                    \
    offset[0] = grid->nx * mpi_rank_x;                                              \
    offset[1] = grid->ny * mpi_rank_y;                                              \
    offset[2] = grid->nz * mpi_rank_z;                                              \
                                                                                    \
    filespace = H5Screate_simple(3, global_size, NULL);                             \
    memspace = H5Screate_simple(3, local_size, NULL);                               \
    dataspace_id;                                                                   \
}
#endif

// read time-averaged fields in HDF5
#ifndef READ_TIME_AVERAGE_FIELDS_HDF5
#define READ_TIME_AVERAGE_FIELDS_HDF5(fbase, DATA_STRUCT, time_step)               \
{                                                                                  \
    hid_t plist_id, file_id, group_id, dset_id;                                    \
    hid_t filespace, memspace, dataspace_id;                                       \
    hsize_t temp_buf_index;                                                        \
    hsize_t global_size[3], local_size[3], offset[3];                              \
    float *temp_buf = (float *)malloc(sizeof(float) *                              \
        (grid->nx) * (grid->ny) * (grid->nz));                                     \
    double el2 = uptime();                                                         \
    OPEN_HDF5_FILE(fbase, "fields", time_step)                                     \
    if (field_dump_flag.ex)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "ex", ex, H5T_NATIVE_FLOAT);                   \
    if (field_dump_flag.ey)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "ey", ey, H5T_NATIVE_FLOAT);                   \
    if (field_dump_flag.ez)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "ez", ez, H5T_NATIVE_FLOAT);                   \
                                                                                   \
    if (field_dump_flag.cbx)                                                       \
        READ_FROM_HDF5(DATA_STRUCT, "cbx", cbx, H5T_NATIVE_FLOAT);                 \
    if (field_dump_flag.cby)                                                       \
        READ_FROM_HDF5(DATA_STRUCT, "cby", cby, H5T_NATIVE_FLOAT);                 \
    if (field_dump_flag.cbz)                                                       \
        READ_FROM_HDF5(DATA_STRUCT, "cbz", cbz, H5T_NATIVE_FLOAT);                 \
                                                                                   \
    el2 = uptime() - el2;                                                          \
    sim_log("TimeHDF5Read: " << el2 << " s");                                      \
                                                                                   \
    double el3 = uptime();                                                         \
                                                                                   \
    free(temp_buf);                                                                \
    H5Sclose(filespace);                                                           \
    H5Sclose(memspace);                                                            \
    H5Pclose(plist_id);                                                            \
    H5Gclose(group_id);                                                            \
    H5Fclose(file_id);                                                             \
                                                                                   \
    el3 = uptime() - el3;                                                          \
    sim_log("TimeHDF5Close: " << el3 << " s");                                     \
}
#endif

// read time-averaged hydro in HDF5
#ifndef READ_TIME_AVERAGE_HYDRO_HDF5
#define READ_TIME_AVERAGE_HYDRO_HDF5(speciesname, fbase, DATA_STRUCT, time_step)    \
{                                                                                   \
    species_t *sp = find_species_name(speciesname, species_list);                   \
    if (!sp)                                                                        \
      ERROR(("Invalid species name: %s", speciesname));                             \
                                                                                    \
    char hname[32];                                                                 \
    sprintf(hname, "hydro_%s", speciesname);                                        \
    hid_t plist_id, file_id, group_id, dset_id;                                     \
    hid_t filespace, memspace, dataspace_id;                                        \
    hsize_t temp_buf_index;                                                         \
    hsize_t global_size[3], local_size[3], offset[3];                               \
    float *temp_buf = (float *)malloc(sizeof(float) *                               \
        (grid->nx) * (grid->ny) * (grid->nz));                                      \
    double el2 = uptime();                                                          \
    OPEN_HDF5_FILE(fbase, hname, time_step)                                         \
                                                                                    \
    if (hydro_dump_flag.jx)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "jx", jx, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.jy)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "jy", jy, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.jz)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "jz", jz, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.rho)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "rho", rho, H5T_NATIVE_FLOAT);                  \
                                                                                    \
    if (hydro_dump_flag.px)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "px", px, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.py)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "py", py, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.pz)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "pz", pz, H5T_NATIVE_FLOAT);                    \
    if (hydro_dump_flag.ke)                                                         \
        READ_FROM_HDF5(DATA_STRUCT, "ke", ke, H5T_NATIVE_FLOAT);                    \
                                                                                    \
    if (hydro_dump_flag.txx)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "txx", txx, H5T_NATIVE_FLOAT);                  \
    if (hydro_dump_flag.tyy)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "tyy", tyy, H5T_NATIVE_FLOAT);                  \
    if (hydro_dump_flag.tzz)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "tzz", tzz, H5T_NATIVE_FLOAT);                  \
                                                                                    \
    if (hydro_dump_flag.tyz)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "tyz", tyz, H5T_NATIVE_FLOAT);                  \
    if (hydro_dump_flag.tzx)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "tzx", tzx, H5T_NATIVE_FLOAT);                  \
    if (hydro_dump_flag.txy)                                                        \
        READ_FROM_HDF5(DATA_STRUCT, "txy", txy, H5T_NATIVE_FLOAT);                  \
                                                                                    \
    el2 = uptime() - el2;                                                           \
    sim_log("TimeHDF5Read: " << el2 << " s");                                       \
                                                                                    \
    double el3 = uptime();                                                          \
                                                                                    \
    free(temp_buf);                                                                 \
    H5Sclose(filespace);                                                            \
    H5Sclose(memspace);                                                             \
    H5Pclose(plist_id);                                                             \
    H5Gclose(group_id);                                                             \
    H5Fclose(file_id);                                                              \
                                                                                    \
    el3 = uptime() - el3;                                                           \
    sim_log("TimeHDF5Close: " << el3 << " s");                                      \
}
#endif

#endif // __TIME_AVERAGE_MASTER_HH__
