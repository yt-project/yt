import yt
import yt.units as u

ds = yt.load('HiresIsolatedGalaxy/DD0044/DD0044')

center = [0.53, 0.53, 0.53]
normal = [0,0,1]
radius = 40*u.kpc
height = 5*u.kpc

disk = ds.disk(center, [0,0,1], radius, height)

profile = yt.create_profile(
    data_source=disk,
    bin_fields=["radius"],
    fields=["cylindrical_tangential_velocity_absolute"],
    n_bins=256,
    units=dict(radius="kpc",
               cylindrical_tangential_velocity_absolute="km/s"),
    logs=dict(radius=False),
    weight_field='cell_mass',
    extrema=dict(radius=(0,40)),
    )

plot = yt.ProfilePlot.from_profiles(profile)

plot.set_log('cylindrical_tangential_velocity_absolute', False)
plot.set_ylim('cylindrical_tangential_velocity_absolute', 60, 160)

plot.save()
