import matplotlib.pyplot as plt
import numpy as np

import yt

# Enable parallelism in the script (assuming it was called with
# `mpirun -np <n_procs>` )
yt.enable_parallelism()

# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.
# The "entropy" field that we will use below depends on the electron number
# density, which is not in these datasets by default, so we assume full
# ionization using the "default_species_fields" kwarg.
ts = yt.load(
    "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0*",
    default_species_fields="ionized",
)

storage = {}

# By using the piter() function, we can iterate on every dataset in
# the TimeSeries object.  By using the storage keyword, we can populate
# a dictionary where the dataset is the key, and sto.result is the value
# for later use when the loop is complete.

# The serial equivalent of piter() here is just "for ds in ts:" .

for store, ds in ts.piter(storage=storage):

    # Create a sphere of radius 100 kpc at the center of the dataset volume
    sphere = ds.sphere("c", (100.0, "kpc"))
    # Calculate the entropy within that sphere
    entr = sphere[("gas", "entropy")].sum()
    # Store the current time and sphere entropy for this dataset in our
    # storage dictionary as a tuple
    store.result = (ds.current_time.in_units("Gyr"), entr)

# Convert the storage dictionary values to a Nx2 array, so the can be easily
# plotted
arr = np.array(list(storage.values()))

# Plot up the results: time versus entropy
plt.semilogy(arr[:, 0], arr[:, 1], "r-")
plt.xlabel("Time (Gyr)")
plt.ylabel("Entropy (ergs/K)")
plt.savefig("time_versus_entropy.png")
