import yt
import yt.units as u

ds = yt.load("HiresIsolatedGalaxy/DD0044/DD0044")

center = [0.53, 0.53, 0.53]
normal = [0, 0, 1]
radius = 40 * u.kpc
height = 2 * u.kpc

disk = ds.disk(center, [0, 0, 1], radius, height)

profile = yt.create_profile(
    data_source=disk,
    bin_fields=[("index", "radius"), ("gas", "velocity_cylindrical_theta")],
    fields=[("gas", "mass")],
    n_bins=256,
    units=dict(radius="kpc", velocity_cylindrical_theta="km/s", mass="Msun"),
    logs=dict(radius=False, velocity_cylindrical_theta=False),
    weight_field=None,
    extrema=dict(radius=(0, 40), velocity_cylindrical_theta=(-250, 250)),
)

plot = yt.PhasePlot.from_profile(profile)
plot.set_cmap(("gas", "mass"), "YlOrRd")

plot.save()
