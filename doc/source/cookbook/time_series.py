from yt.mods import *
import glob
import matplotlib.pyplot as plt

keV = 1.16044e7
mue = 1.0/0.875
m_p = 1.673e-24
mtt = -2./3.

# Glob for a list of filenames, then sort them
fns = glob.glob("GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0[0-6][0-9]0")
fns.sort()

# Construct the time series object

ts = DatasetSeries.from_filenames(fns)

storage = {}

# We use the piter() method here so that this can be run in parallel.
# Alternately, you could just iterate "for ds in ts:" and directly append to
# times and entrs.
for sto, ds in ts.piter(storage=storage):
    sphere = ds.sphere("c", (100., "kpc"))
    temp = sphere["temperature"]/keV
    dens = sphere["density"]/(m_p*mue)
    mgas = sphere["cell_mass"]
    entr = (temp*(dens**mtt)*mgas).sum()/mgas.sum() 
    sto.result = (ds.current_time, entr)

times = []
entrs = []
for k in storage:
    t, e = storage[k]
    times.append(t)
    entrs.append(e)

if is_root():
    plt.semilogy(times, entrs, 'x-')
    plt.xlabel("Time")
    plt.ylabel("Entropy")
    plt.savefig("time_versus_entropy.png")
