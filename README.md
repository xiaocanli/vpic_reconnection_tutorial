# Running Reconnection Simulations Using VPIC
The tutorial is based on magnetic reconnection simulations. The simulation decks are in [decks](decks). We will assume you can access the NERSC's Cori or TACC's Frontera, but the procedure can be generalized to other machines without too much additional effort.

The old version of VPIC uses binary (N-to-N) output (it is still the default for the master branch). The newer versions now support HDF5 file output. We will recommend using the HDF5 version of VPIC.

Check out how to compile VPIC HDF5 branch and do a test run at [HDF5_version_of_VPIC](HDF5_version_of_VPIC.md). After that, you should be able to produce some nice pictures of the reconnection layer.

If you decide to modify the reconnection deck, see how to do that at [modify_reconnection_deck](modify_reconnection_deck.md). Make sure to recompile the deck after you change some variables.

If you hope to tracer some particles in high cadence to plot some particle trajectories as the simulation is running, checkout [particle_tracking](particle_tracking.md).

Check out how to dump particle data at [modify_reconnection_deck](modify_reconnection_deck.md). You can plot some phase space diagram using the particle data. Check out a simple script [phase_diagram.py](decks/mime1836_sigmaic256_bg00/phase_diagram.py) to do that.
