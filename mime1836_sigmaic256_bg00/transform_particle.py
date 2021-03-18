"""
Transform particles from local frame to global frame
"""
import argparse
import errno
import math
import multiprocessing
import os

import h5py
import numpy as np
from joblib import Parallel, delayed

# give some PIC simulation parameters here
lx_pic = 250.0  # in electron inertial length
ly_pic = 25.0
lz_pic = 125.0
nx_pic, ny_pic, nz_pic = 1024, 1, 512
topo_x, topo_y, topo_z = 32, 1, 2
particle_interval = 517

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_var(group, dset_name, dset_size):
    """Read data from a HDF5 group

    Args:
        group: HDF5 group
        dset_name: the dataset name
        dset_size: the size of the data
    """
    dset = group[dset_name]
    fdata = np.zeros(dset_size, dtype=dset.dtype)
    dset.read_direct(fdata)
    return fdata


def frame_transformation(fname, meta_name, tindex, particle_name):
    """Transform particles from local frame to global frame

    Args:
        fname: the name of the H5Part file
        meta_name: the meta filename
        tindex: time index
        particle_name: name of the particle
    """
    with h5py.File(fname, 'r') as fh:
        group_name = 'Step#' + str(tindex)
        group = fh[group_name]
        dset = group['dX']
        nptl, = dset.shape
        ptl = {}
        for dset in group:
            dset = str(dset)
            ptl[str(dset)] = read_var(group, dset, nptl)

    with h5py.File(meta_name, 'r') as fh:
        group_name = 'Step#' + str(tindex)
        group = fh[group_name]
        dset = group['np_local']
        mpi_size, = dset.shape
        meta_data = {}
        for dset in group:
            dset = str(dset)
            meta_data[str(dset)] = read_var(group, dset, mpi_size)

    # grid sizes
    dx_pic = lx_pic / nx_pic
    dy_pic = ly_pic / ny_pic
    dz_pic = lz_pic / nz_pic

    # local domain dimensions
    nx_mpi = nx_pic // topo_x + 2
    ny_mpi = ny_pic // topo_y + 2
    nz_mpi = nz_pic // topo_z + 2

    zcell = ptl["i"] // (nx_mpi * ny_mpi)  # [1, nx_mpi - 1]
    ycell = (ptl["i"] % (nx_mpi * ny_mpi)) // nx_mpi
    xcell = ptl["i"] % nx_mpi

    ptl["X"] = np.float32(0.5 * (ptl["dX"] + 1.0) + (xcell - 1) * dx_pic)
    ptl["Y"] = np.float32(0.5 * (ptl["dY"] + 1.0) + (ycell - 1) * dy_pic)
    ptl["Z"] = np.float32(0.5 * (ptl["dZ"] + 1.0) + (zcell - 1) * dz_pic)

    ptl.pop("dX", None)
    ptl.pop("dY", None)
    ptl.pop("dZ", None)

    np_local = meta_data["np_local"]
    np_offset = np.zeros(mpi_size + 1, dtype=np.int64)
    np_offset[1:] = np.cumsum(np_local)
    for rank in range(mpi_size):
        start = np_offset[rank]
        end = np_offset[rank+1]
        ptl["X"][start:end] += meta_data["x0"][rank]
        ptl["Y"][start:end] += meta_data["y0"][rank]
        ptl["Z"][start:end] += meta_data["z0"][rank]

    fdir = "particle/particle_transformed/"
    mkdir_p(fdir)
    # fname = fdir + particle_name + "_" + str(tindex) + ".h5part"
    # with h5py.File(fname, 'w') as fh:
    #     group_name = "Step#" + str(tindex)
    # fname = fdir + particle_name + ".h5part"
    fname = fdir + particle_name + "_reduced.h5part"
    with h5py.File(fname, 'a') as fh:
        group_name = "Step#" + str(tindex//particle_interval - 1)
        group = fh.create_group(group_name)
        for dset in ptl.keys():
            group.create_dataset(dset, (nptl, ), data=ptl[dset])


def get_cmd_args():
    """Get command line arguments """
    default_run_name = "test_2d_1"
    default_run_dir = ("/net/scratch3/xiaocanli/reconnection/Cori_runs/" +
                       default_run_name + "/")
    parser = argparse.ArgumentParser(description="Transform particle")
    parser.add_argument("--run_dir", action="store", default=default_run_dir,
                        help="PIC run directory")
    parser.add_argument("--run_name", action="store", default=default_run_name,
                        help="PIC run name")
    parser.add_argument("--output_dir", action="store", default=default_run_dir,
                        help="Output root directory for transformed particles")
    parser.add_argument('--tframe', action="store", default='0', type=int,
                        help='Time frame')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='Starting time frame')
    parser.add_argument('--tend', action="store", default='8', type=int,
                        help='Ending time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to analyze multiple frames using a time loop')
    parser.add_argument("--species", action="store", default="electron",
                        help="particle species")
    return parser.parse_args()


def analysis_single_frame(args):
    """Analysis for single time frame
    """
    print("Time frame: %d" % args.tframe)
    fname = (args.run_dir + "particle//T." + str(args.tframe) +
             "/" + args.species + "_" + str(args.tframe) + ".h5part")
    meta_name = (args.run_dir + "particle/T." + str(args.tframe) +
                 "/grid_metadata_" + args.species + "_" +
                 str(args.tframe) + ".h5part")
    frame_transformation(fname, meta_name, args.tframe, args.species)


def process_input(args, tframe):
    """process one time frame"""
    print("Time frame: %d" % tframe)
    fname = (args.run_dir + "particle//T." + str(tframe) +
             "/" + args.species + "_" + str(args.tframe) + ".h5part")
    meta_name = (args.run_dir + "particle/T." + str(tframe) +
                 "/grid_metadata_" + args.species + "_" +
                 str(args.tframe) + ".h5part")
    frame_transformation(fname, meta_name, tframe, args.species)


def analysis_multi_frames(args):
    """Analysis for multiple time frames
    """
    tframes = range(args.tstart, args.tend + 1, particle_interval)
    if args.time_loop:
        for tframe in tframes:
            args.tframe = tframe
            analysis_single_frame(args)
    else:
        ncores = multiprocessing.cpu_count()
        # ncores = 8
        Parallel(n_jobs=ncores)(delayed(process_input)(args, tframe)
                                for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    if args.multi_frames:
        analysis_multi_frames(args)
    else:
        analysis_single_frame(args)


if __name__ == "__main__":
    main()
