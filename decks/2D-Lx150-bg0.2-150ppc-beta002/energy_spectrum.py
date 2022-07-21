"""
#!/usr/bin/env python3
"""
from __future__ import print_function

import argparse
import itertools
import json
import multiprocessing

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# from joblib import Parallel, delayed
# from matplotlib.colors import LogNorm

# mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# mpl.rc('text', usetex=True)
# mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

def read_spectrum(species, tframe, nbins):
    """Read energy spectrum for a HDF5 file

    Args:
        species: particle species, "electron", "ion" or others
        tframe: time frame
        nbins: number of energy bins
    """
    fname = ("spectrum/T." + str(tframe) + "/spectrum_" +
             species + "_" + str(tframe) + ".h5part")
    print(fname)
    with h5py.File(fname, 'r') as fh:
        # group_name = "Step#" + str(tframe)
        # group = fh[group_name]
        dset = fh["spectrum"]
        sz, = dset.shape
        fdata = np.zeros(sz, dtype=dset.dtype)
        dset.read_direct(fdata)

    ntot = nbins + 3  # including bx, by, and bz
    ebins = np.logspace(-5, 5, nbins)
    print(fdata.shape[0])
    # fdata = fdata.reshape((fdata.shape[0]//ntot, ntot))
    # spectrum = np.sum(fdata[:, 3:], axis=0) / np.gradient(ebins)
    # plt.loglog(ebins, spectrum)
    # plt.show()

def plot_spectrum(species):
    fig = plt.figure(figsize=[7, 8])
    w1, h1 = 0.83, 0.32
    xs, ys = 0.96 - w1, 0.80 - h1
    ax = fig.add_axes([xs, ys, w1, h1])
    # for tframe in range(0, 88680, 2217):
    for tframe in range(88680, 88681, 2217):
        fname = "spectrum_combined/spectrum_" + species + "_" + str(tframe) + ".dat"
        spect = np.fromfile(fname, dtype=np.float32)
        ndata, = spect.shape
        ebins = np.logspace(-6, 4, ndata-3)
        spect[3:] /= np.gradient(ebins)
        ax.loglog(ebins, spect[3:])
    fpower = 1E9*ebins**-4
    ax.loglog(ebins, fpower)
    ax.set_xlim([1E-3, 2E1])
    ax.set_ylim([1E3, 1E15])
    plt.show()


def get_cmd_args():
    """Get command line arguments
    """
    default_pic_run = 'bb_test_integrated'
    default_pic_run_dir = ('/global/cscratch1/sd/xiaocan/Cori_tests' +
                           default_pic_run + '/')
    parser = argparse.ArgumentParser(description='Particle energy spectrum')
    parser.add_argument('--pic_run', action="store",
                        default=default_pic_run, help='PIC run name')
    parser.add_argument('--pic_run_dir', action="store",
                        default=default_pic_run_dir, help='PIC run directory')
    parser.add_argument('--species', action="store",
                        default="electron", help='Particle species')
    parser.add_argument('--tframe', action="store", default='0', type=int,
                        help='Time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    return parser.parse_args()


def analysis_single_frames(plot_config, args):
    """Analysis for multiple time frames
    """
    # read_spectrum(args.species, args.tframe, 1000)
    plot_spectrum('e')
    plot_spectrum('i')


def process_input(plot_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    pass


def analysis_multi_frames(plot_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    ncores = multiprocessing.cpu_count()
    # Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, args, tframe)
    #                         for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    plot_config = {}
    if args.multi_frames:
        analysis_multi_frames(plot_config, args)
    else:
        analysis_single_frames(plot_config, args)


if __name__ == "__main__":
    main()
