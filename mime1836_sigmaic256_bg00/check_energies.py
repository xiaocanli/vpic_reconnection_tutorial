import numpy as np
import matplotlib.pylab as plt

fdata = np.genfromtxt("rundata/energies", skip_header=3)
ene_electric = np.sum(fdata[:, 1:4], axis=1)
ene_magnetic = np.sum(fdata[:, 4:7], axis=1)
ene_ion = fdata[:, 7]
ene_electron = fdata[:, 8]
ene_tot = np.sum(fdata[:, 1:], axis=1)

print("Total energy: %f" % (ene_tot[-1]/ene_tot[0]))

enorm = ene_magnetic[0]
# plt.plot(ene_electric/enorm)
# plt.plot(ene_magnetic/enorm)
# plt.plot(ene_ion/enorm)
# plt.plot(ene_electron/enorm)

# plt.plot(ene_tot / ene_tot[0])

# fdata = np.genfromtxt("energies_cfl077_nppc100", skip_header=3)
# ene_tot = np.sum(fdata[:, 1:], axis=1)
# plt.plot(ene_tot / ene_tot[0], label="cfl077-nppc100")

# fdata = np.genfromtxt("energies_cfl070_nppc100", skip_header=3)
# ene_tot = np.sum(fdata[:, 1:], axis=1)
# plt.plot(ene_tot / ene_tot[0], label="cfl070-nppc100")

# fdata = np.genfromtxt("energies_cfl070_nppc200", skip_header=3)
# ene_tot = np.sum(fdata[:, 1:], axis=1)
# plt.plot(ene_tot / ene_tot[0], label="cfl070-nppc200")

# plt.grid()

# plt.legend()

# plt.show()
