def plot_trajectory(plot_config, show_plot=True):
    """Plot particle trajectory
    """
    species = plot_config["species"]
    pic_run = plot_config["pic_run"]
    tracer_dir = "/net/scratch3/xiaocanli/vpic-sorter/data/trans-relativistic/"
    picinfo_fname = '../data/pic_info/pic_info_' + pic_run + '.json'
    pic_info = read_data_from_json(picinfo_fname)
    dtwpe_tracer = pic_info.dtwpe * pic_info.tracer_interval
    fname = tracer_dir + pic_run + "/electrons_ntraj500_1emax.h5p"
    fh = h5py.File(fname, 'r')
    particle_tags = list(fh.keys())
    nptl = len(particle_tags)
    print("Total number of particles: %d" % nptl)
    group = fh[particle_tags[0]]
    dset = group['dX']
    nframes, = dset.shape
    ptl = {}
    for dset_name in group:
        dset = group[dset_name]
        ptl[str(dset_name)] = np.zeros(dset.shape, dset.dtype)
        dset.read_direct(ptl[str(dset_name)])
    gamma0 = 1.0 / np.sqrt(1 + ptl["Ux"]**2 + ptl["Uy"]**2 + ptl["Uz"]**2)
    ttracer = np.arange(0, nframes) * dtwpe_tracer
    tmin, tmax = ttracer[0], ttracer[-1]

    semilogy = False

    if semilogy:
        img_dir = '../img/trans_relativistic/tracer_traj_log/' + pic_run + '/'
    else:
        img_dir = '../img/trans_relativistic/tracer_traj/' + pic_run + '/'
    mkdir_p(img_dir)
    if species in ["e", "electron"]:
        sname = "electron"
    else:
        sname = "Ion"

    for iptl in range(nptl):
        # for iptl in range(1, 2):
        print("Particle: %d of %d" % (iptl, nptl))
        group = fh[particle_tags[iptl]]
        for dset_name in group:
            dset = group[dset_name]
            dset.read_direct(ptl[str(dset_name)])
        gamma = np.sqrt(1 + ptl["Ux"]**2 + ptl["Uy"]**2 + ptl["Uz"]**2)
        dgamma = gamma - gamma0
        igamma = 1.0 / gamma
        vx = ptl["Ux"] * igamma
        vy = ptl["Uy"] * igamma
        vz = ptl["Uz"] * igamma
        x = ptl["dX"]
        y = ptl["dY"]
        z = ptl["dZ"]
        ex = ptl["Ex"]
        ey = ptl["Ey"]
        ez = ptl["Ez"]
        bx = ptl["Bx"]
        by = ptl["By"]
        bz = ptl["Bz"]
        edotb = ex * bx + ey * by + ez * bz
        ib2 = 1.0 / (bx**2 + by**2 + bz**2)
        eparax = edotb * bx * ib2
        eparay = edotb * by * ib2
        eparaz = edotb * bz * ib2
        wtot = np.cumsum(-(ex * vx + ey * vy + ez * vz)) * dtwpe_tracer
        wpara = np.cumsum(-(eparax * vx + eparay * vy +
                            eparaz * vz)) * dtwpe_tracer
        wperp = wtot - wpara
        fig = plt.figure(figsize=[5, 3.5])
        rect = [0.14, 0.16, 0.82, 0.8]
        ax = fig.add_axes(rect)
        COLORS = palettable.tableau.Tableau_10.mpl_colors
        ax.set_prop_cycle('color', COLORS)
        if semilogy:
            ax.semilogy(ttracer, wpara, linewidth=2, label=r'$W_\parallel$')
            ax.semilogy(ttracer, wperp, linewidth=2, label=r'$W_\perp$')
            ax.semilogy(ttracer, dgamma, linewidth=2, label=r'$\Delta\gamma$')
        else:
            ax.plot(ttracer, wpara, linewidth=2, label=r'$W_\parallel$')
            ax.plot(ttracer, wperp, linewidth=2, label=r'$W_\perp$')
            ax.plot(ttracer, dgamma, linewidth=2, label=r'$\Delta\gamma$')
        ax.set_xlim([tmin, tmax])
        ax.tick_params(bottom=True, top=True, left=True, right=True)
        ax.tick_params(axis='x', which='minor', direction='in')
        ax.tick_params(axis='x', which='major', direction='in')
        ax.tick_params(axis='y', which='minor', direction='in')
        ax.tick_params(axis='y', which='major', direction='in')
        ax.set_xlabel(r'$t\omega_{pe}$', fontsize=16)
        ax.set_ylabel('Energy change', fontsize=16)
        ax.tick_params(labelsize=12)
        ax.set_xlim([0, 1.5E4])
        if semilogy:
            ax.set_ylim([1E-1, 3E3])
            ax.legend(loc=4,
                      prop={'size': 12},
                      ncol=1,
                      shadow=False,
                      fancybox=False,
                      frameon=False)
        else:
            ax.legend(loc=6,
                      prop={'size': 12},
                      ncol=1,
                      shadow=False,
                      fancybox=False,
                      frameon=False)
        fname = img_dir + sname + "_tracer_" + str(iptl) + ".pdf"
        fig.savefig(fname)

        plt.close()
        # plt.show()

    fh.close()
