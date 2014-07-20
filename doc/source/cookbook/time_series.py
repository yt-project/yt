import yt
import glob
import matplotlib.pyplot as plt

# Glob for a list of filenames, then sort them
fns = glob.glob("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0*")
fns.sort()

# Construct the time series object
ts = yt.DatasetSeries.from_filenames(fns)

storage = {}

# We use the piter() method here so that this can be run in parallel.
# Alternately, you could just iterate "for ds in ts:" and directly append to
# times and entrs.
for sto, ds in ts.piter(storage=storage):
    sphere = ds.sphere("c", (100., "kpc"))
    entr = sphere["entropy"].sum()
    sto.result = (ds.current_time.in_units('Gyr'), entr)


# Store these values in a couple of lists
times = []
entrs = []
for k in storage:
    t, e = storage[k]
    times.append(t)
    entrs.append(e)


# Plot up the results

plt.semilogy(times, entrs, '-')
plt.xlabel("Time (Gyr)")
plt.ylabel("Entropy (ergs/K)")
plt.savefig("time_versus_entropy.png")
