// This requires parallel HDF5.
// Compile on Cori: CC -o compress_hdf5 compress_hdf5.cc
// Before running the code, create data_smooth, and set stripe if using Lustre
// file system. For example, "stripe_medium data_smooth" on Cori.
//
#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>
#include <string.h>
#include "hdf5.h"

// Variables that might need to be modified
const int nx = 3072; // Change to your simulation sizes
const int ny = 1536;
const int nz = 1280;
const int reduce = 2;
const int topox = 32; // MPI topology
const int topoy = 1;
const int topoz = 1;
const int tinterval = 2217;
const int nframes = 41;
int nx_l, ny_l, nz_l, nxr_l, nyr_l, nzr_l;
void reduce_data(std::string filename, std::string groupname,
    std::vector<std::string> var_names, std::string output_filename,
    float *input, float *output, hsize_t subsizes[3], hsize_t starts[3],
    hsize_t sizes_out[3], hsize_t subsizes_out[3], hsize_t starts_out[3]);

int main(int argc, char **argv)
{
  int mpi_rank, mpi_size;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  int nxr = nx / reduce;
  int nyr = ny / reduce;
  int nzr = nz / reduce;

  double t0 = MPI_Wtime();

  if ((topox * topoy * topoz) != mpi_size ||
      (nx / topox / reduce) * topox * reduce != nx ||
      (ny / topoy / reduce) * topoy * reduce != ny ||
      (nz / topoz / reduce) * topoz * reduce != nz) {
    if (mpi_rank == 0) printf("ERROR: %s\n", "Wrong MPI Topology");
    MPI_Finalize();
    exit(1);
  }

  int mpi_rankx = mpi_rank / (topoy * topoz);
  int mpi_ranky = (mpi_rank % (topoy * topoz)) / topoz;
  int mpi_rankz = mpi_rank % topoz;

  // input sizes
  nx_l = nx / topox;
  ny_l = ny / topoy;
  nz_l = nz / topoz;
  hsize_t sizes[3] = {nx, ny, nz};
  hsize_t subsizes[3] = {nx_l, ny_l, nz_l};
  hsize_t starts[3] = {nx_l*mpi_rankx, ny_l*mpi_ranky, nz_l*mpi_rankz};
  // output sizes
  nxr_l = nxr / topox;
  nyr_l = nyr / topoy;
  nzr_l = nzr / topoz;
  hsize_t sizes_out[3] = {nxr, nyr, nzr};
  hsize_t subsizes_out[3] = {nxr_l, nyr_l, nzr_l};
  hsize_t starts_out[3] = {nxr_l*mpi_rankx, nyr_l*mpi_ranky, nzr_l*mpi_rankz};

  // input and output data
  float *input = (float *)malloc(nx_l * ny_l * nz_l * sizeof(float));
  float *output = (float *)malloc(nxr_l * nyr_l * nzr_l * sizeof(float));

  static std::vector<std::string> field_names = {"cbx", "cby", "cbz", "ex", "ey", "ez"};
  static std::vector<std::string> hydro_names = {"jx", "jy", "jz", "ke", "px", "py", "pz",
                                                 "rho", "txx", "txy", "tyy", "tyz", "tzx", "tzz"};

  for (int iframe = 20; iframe != 21; ++iframe) {
    if (mpi_rank == 0) std::cout << iframe << " of " << nframes << std::endl;
    double tpre = MPI_Wtime();
    int tindex = iframe * tinterval;
    // EMF
    std::string filename = std::string("field_hdf5/T.") +
      std::to_string(tindex) + "/fields_" + std::to_string(tindex) + ".h5";
    std::string groupname = std::string("Timestep_") + std::to_string(tindex);
    std::string output_filename = std::string("data_smooth/fields_") + std::to_string(tindex) + ".h5";
    reduce_data(filename, groupname, field_names, output_filename, input, output,
        subsizes, starts, sizes_out, subsizes_out, starts_out);
    // Electron hydro
    filename = std::string("hydro_hdf5/T.") + std::to_string(tindex) +
      "/hydro_electron_" + std::to_string(tindex) + ".h5";
    output_filename = std::string("data_smooth/hydro_electron_") + std::to_string(tindex) + ".h5";
    reduce_data(filename, groupname, hydro_names, output_filename, input, output,
        subsizes, starts, sizes_out, subsizes_out, starts_out);
    // Ion hydro
    filename = std::string("hydro_hdf5/T.") + std::to_string(tindex) +
      "/hydro_ion_" + std::to_string(tindex) + ".h5";
    output_filename = std::string("data_smooth/hydro_ion_") + std::to_string(tindex) + ".h5";
    reduce_data(filename, groupname, hydro_names, output_filename, input, output,
        subsizes, starts, sizes_out, subsizes_out, starts_out);
    double tpos = MPI_Wtime();
    if(mpi_rank == 0)
      printf("Time for one step: [%f]s \n", (tpos - tpre));
  }

  free(input);
  free(output);

  double t1 = MPI_Wtime();
  if(mpi_rank == 0) {
    printf("Overall time is [%f]s \n", (t1 - t0));
  }

  MPI_Finalize();

  return 0;
}

/* ----------------------------------------------------------------------
 * read data from file, reduce the data, and save the data
 * ---------------------------------------------------------------------- */
void reduce_data(std::string filename, std::string groupname,
    std::vector<std::string> var_names, std::string output_filename,
    float *input, float *output, hsize_t subsizes[3], hsize_t starts[3],
    hsize_t sizes_out[3], hsize_t subsizes_out[3], hsize_t starts_out[3])
{
  MPI_Info info;
  MPI_Info_create(&info);
  MPI_Info_set(info, "romio_cb_read", "enable");
  MPI_Info_set(info, "romio_cb_write", "enable");

  // open input and output file and group
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
  hid_t output_file_id = H5Fcreate(output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  hid_t group_id = H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT);
  hid_t output_group_id = H5Gcreate(output_file_id, groupname.c_str(),
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t memspace = H5Screate_simple(3, subsizes, NULL);
  hid_t output_filespace = H5Screate_simple(3, sizes_out, NULL);
  hid_t output_memspace = H5Screate_simple(3, subsizes_out, NULL);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  for (auto var : var_names) {
    // read
    hid_t dset_id = H5Dopen(group_id, var.c_str(), H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, subsizes, NULL);
    H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, input);
    H5Sclose(filespace);
    H5Dclose(dset_id);
    // reduce
    int ir, jr, kr;
    float idv = 1.0 / (reduce*reduce*reduce);
    memset(output, 0, nxr_l * nyr_l * nzr_l * sizeof(float));
    for (int i = 0; i != nx_l; ++i) {
      ir = i / reduce;
      for (int j = 0; j != ny_l; ++j) {
        jr = j / reduce;
        for (int k = 0; k != nz_l; ++k) {
          kr = k / reduce;
          output[kr+nzr_l*(jr+ir*nyr_l)] += input[k+nz_l*(j+i*ny_l)] * idv;
        }
      }
    }
    // write
    dset_id = H5Dcreate(output_group_id, var.c_str(), H5T_NATIVE_FLOAT,
        output_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(output_filespace, H5S_SELECT_SET, starts_out, NULL, subsizes_out, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, output_memspace, output_filespace, plist_id, output);
    H5Dclose(dset_id);
  }

  MPI_Info_free(&info);

  H5Pclose(plist_id);
  H5Sclose(output_memspace);
  H5Sclose(output_filespace);
  H5Gclose(output_group_id);
  H5Fclose(output_file_id);
  H5Sclose(memspace);
  H5Gclose(group_id);
  H5Fclose(file_id);
}
