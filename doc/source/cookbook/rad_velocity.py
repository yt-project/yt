import matplotlib.pyplot as plt

import yt

ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

# Get the first sphere
sp0 = ds.sphere(ds.domain_center, (500.0, "kpc"))

# Compute the bulk velocity from the cells in this sphere
bulk_vel = sp0.quantities.bulk_velocity()


# Get the second sphere
sp1 = ds.sphere(ds.domain_center, (500.0, "kpc"))

# Set the bulk velocity field parameter
sp1.set_field_parameter("bulk_velocity", bulk_vel)

# Radial profile without correction

rp0 = yt.create_profile(
    sp0,
    ("index", "radius"),
    ("gas", "radial_velocity"),
    units={("index", "radius"): "kpc"},
    logs={("index", "radius"): False},
)

# Radial profile with correction for bulk velocity

rp1 = yt.create_profile(
    sp1,
    ("index", "radius"),
    ("gas", "radial_velocity"),
    units={("index", "radius"): "kpc"},
    logs={("index", "radius"): False},
)

# Make a plot using matplotlib

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(
    rp0.x.value,
    rp0[("gas", "radial_velocity")].in_units("km/s").value,
    rp1.x.value,
    rp1[("gas", "radial_velocity")].in_units("km/s").value,
)

ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{v_r\ (km/s)}$")
ax.legend(["Without Correction", "With Correction"])

fig.savefig(f"{ds}_profiles.png")
