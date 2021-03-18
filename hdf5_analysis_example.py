import math

import h5py
import matplotlib.pylab as plt
import numpy as np
from matplotlib import rc

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

plt.rcParams['figure.dpi'] = 100

tableau_colors=['tab:blue', 'tab:orange', 'tab:green',
                'tab:red', 'tab:purple', 'tab:brown',
                'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


pic_run = "driving_test_dipole"
pic_run_dir = "/global/cscratch1/sd/xiaocan/tail_problem/" + pic_run + "/"
vpic_info = get_vpic_info(pic_run_dir)

def get_vpic_info(pic_run_dir):
    """Get information of the VPIC simulation
    """
    with open(pic_run_dir + '/info') as f:
        content = f.readlines()
    f.close()
    vpic_info = {}
    for line in content[1:]:
        if "=" in line:
            line_splits = line.split("=")
        elif ":" in line:
            line_splits = line.split(":")

        tail = line_splits[1].split("\n")
        vpic_info[line_splits[0].strip()] = float(tail[0])
    return vpic_info


def plot_jy(tframe, yslice=0, show_plot=True):
    """Plot the y-component of the current density
    """
    fields_interval = int(vpic_info["fields_interval"])
    tindex = fields_interval * tframe
    smime = math.sqrt(vpic_info["mi/me"])
    lx_de = vpic_info["Lx/di"] * smime
    lz_de = vpic_info["Lz/di"] * smime
    xmin, xmax = 0, lx_de
    zmin, zmax = -0.5 * lz_de, 0.5 * lz_de

    fname = (pic_run_dir + "hydro_hdf5/T." + str(tindex) +
             "/hydro_electron_" + str(tindex) + ".h5")
    with h5py.File(fname, 'r') as fh:
        group = fh["Timestep_" + str(tindex)]
        dset = group["jy"]
        jey = dset[:, yslice, :]

    fname = (pic_run_dir + "hydro_hdf5/T." + str(tindex) +
             "/hydro_ion_" + str(tindex) + ".h5")
    with h5py.File(fname, 'r') as fh:
        group = fh["Timestep_" + str(tindex)]
        dset = group["jy"]
        jiy = dset[:, yslice, :]

    nx = int(vpic_info["nx"])
    nz = int(vpic_info["nz"])
    xgrid = np.linspace(xmin, xmax, nx)
    zgrid = np.linspace(zmin, zmax, nz)

    len0 = 10
    fig = plt.figure(figsize=[len0, len0*lz_de/lx_de])
    rect = [0.12, 0.14, 0.78, 0.78]
    ax = fig.add_axes(rect)
    jy = np.squeeze(jey + jiy)
    im1 = ax.imshow(jy.T,
                    extent=[xmin, xmax, zmin, zmax],
                    vmin=-0.06, vmax=0.06,
                    cmap=plt.cm.coolwarm, aspect='auto',
                    origin='lower', interpolation='bicubic')
    # Magnetic field lines
    fname = (pic_run_dir + "field_hdf5/T." + str(tindex) +
             "/fields_" + str(tindex) + ".h5")
    bvec = {}
    with h5py.File(fname, 'r') as fh:
        group = fh["Timestep_" + str(tindex)]
        for var in ["cbx", "cbz"]:
            dset = group[var]
            bvec[var] = dset[:, 0, :]
    xmesh, zmesh = np.meshgrid(xgrid, zgrid)
    xmesh_r, zmesh_r = np.meshgrid(xgrid[::16], zgrid[::16])
    start_points = np.vstack([xmesh_r.flatten(), zmesh_r.flatten()]).T
    ax.streamplot(xmesh, zmesh, np.squeeze(bvec["cbx"]).T,
                  np.squeeze(bvec["cbz"]).T, color='k',
                  linewidth=0.5, density=2)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([zmin, zmax])
    ax.set_xlabel(r'$x/d_e$', fontsize=20)
    ax.set_ylabel(r'$z/d_e$', fontsize=20)
    ax.tick_params(labelsize=16)
    rect_cbar = np.copy(rect)
    rect_cbar[0] += rect[2] + 0.01
    rect_cbar[2] = 0.02
    rect_cbar[1] += rect[3] * 0.25
    rect_cbar[3] = rect[3] * 0.5
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(im1, cax=cbar_ax, extend='both')
    cbar_ax.set_title(r'$j_y$', fontsize=20)
    cbar.ax.tick_params(labelsize=12)
    twpe = math.ceil(tindex * vpic_info["dt*wpe"] / 0.1) * 0.1
    text1 = r'$t\omega_{pe}=' + ("{%0.0f}" % twpe) + '$'
    fig.suptitle(text1, fontsize=20)
#     img_dir = '../img/rate_problem/absj/' + pic_run + '/'
#     mkdir_p(img_dir)
#     fname = img_dir + "absj_" + str(tframe) + ".jpg"
#     fig.savefig(fname, dpi=200)


plot_jy(tframe, yslice=0)
