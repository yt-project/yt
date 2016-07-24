import yt

# Create a time-series object.
sim = yt.simulation("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo")
sim.get_time_series(redshifts=[5, 4, 3, 2, 1, 0])

# Lists to hold profiles, labels, and plot specifications.
profiles = []
labels = []
plot_specs = []

# Loop over each dataset in the time-series.
for ds in sim:
    # Create a data container to hold the whole dataset.
    ad = ds.all_data()
    # Create a 1d profile of density vs. temperature.
    profiles.append(yt.create_profile(ad, ["density"],
                                   fields=["temperature"]))
    # Add labels and linestyles.
    labels.append("z = %.2f" % ds.current_redshift)
    plot_specs.append(dict(linewidth=2, alpha=0.7))

# Create the profile plot from the list of profiles.
plot = yt.ProfilePlot.from_profiles(profiles, labels=labels,
                                 plot_specs=plot_specs)
# Save the image.
plot.save()
