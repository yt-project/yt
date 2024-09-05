import numpy as np
from matplotlib import pyplot as plt

import yt
from yt.data_objects.particle_filters import add_particle_filter


def formed_star(pfilter, data):
    filter = data["all", "creation_time"] > 0
    return filter


add_particle_filter(
    "formed_star", function=formed_star, filtered_type="all", requires=["creation_time"]
)

filename = "IsolatedGalaxy/galaxy0030/galaxy0030"

ds = yt.load(filename)
ds.add_particle_filter("formed_star")
ad = ds.all_data()
masses = ad["formed_star", "particle_mass"].in_units("Msun")
formation_time = ad["formed_star", "creation_time"].in_units("yr")

time_range = [0, 5e8]  # years
n_bins = 1000
hist, bins = np.histogram(
    formation_time,
    bins=n_bins,
    range=time_range,
)
inds = np.digitize(formation_time, bins=bins)
time = (bins[:-1] + bins[1:]) / 2

sfr = np.array(
    [masses[inds == j + 1].sum() / (bins[j + 1] - bins[j]) for j in range(len(time))]
)
sfr[sfr == 0] = np.nan

plt.plot(time / 1e6, sfr)
plt.xlabel("Time  [Myr]")
plt.ylabel(r"SFR  [M$_\odot$ yr$^{-1}$]")
plt.savefig("filter_sfr.png")
