import yt
import yt.units as u

ds = yt.load('HiresIsolatedGalaxy/DD0044/DD0044')

center = [0.53, 0.53, 0.53]
normal = [0,0,1]
radius = 40*u.kpc
height = 2*u.kpc

disk = ds.disk(center, [0,0,1], radius, height)

profile = yt.create_profile(
    data_source=disk,
    bin_fields=["radius", "cylindrical_tangential_velocity"],
    fields=["cell_mass"],
    n_bins=256,
    units=dict(radius="kpc",
               cylindrical_tangential_velocity="km/s",
               cell_mass="Msun"),
    logs=dict(radius=False,
              cylindrical_tangential_velocity=False),
    weight_field=None,
    extrema=dict(radius=(0,40),
                 cylindrical_tangential_velocity=(-250, 250)),
    )

plot = yt.PhasePlot.from_profile(profile)
plot.set_cmap("cell_mass", "YlOrRd")

plot.save()
