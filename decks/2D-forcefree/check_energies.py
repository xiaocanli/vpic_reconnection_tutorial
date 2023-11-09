import numpy as np
import matplotlib.pylab as plt

fdata = np.genfromtxt("rundata/energies", skip_header=3)
ene_electric = np.sum(fdata[:, 1:4], axis=1)
ene_magnetic = np.sum(fdata[:, 4:7], axis=1)
ene_ion = fdata[:, 7]
ene_electron = fdata[:, 8]
ene_tot = np.sum(fdata[:, 1:], axis=1)
enorm = ene_magnetic[0]

print("Total energy: %f" % (ene_tot[-1]/ene_tot[0]))

# enorm = ene_magnetic[0]
plt.plot(ene_electric/enorm)
plt.plot(ene_magnetic/enorm)
plt.plot(ene_ion/enorm)
plt.plot(ene_electron/enorm)
plt.savefig("ene.pdf")

# plt.show()
