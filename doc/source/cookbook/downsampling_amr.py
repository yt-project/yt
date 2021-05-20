import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# The maximum refinement level of this dataset is 8
print(ds.max_level)

# If we ask for *all* of the AMR data, we get back field
# values sampled at about 3.6 million AMR zones
ad = ds.all_data()
print(ad["gas", "density"].shape)

# Let's only sample data up to AMR level 2
ad.max_level = 2

# Now we only sample from about 200,000 zones
print(ad["gas", "density"].shape)

# Note that this includes data at level 2 that would
# normally be masked out. There aren't any "holes" in
# the downsampled AMR mesh, the volume still sums to
# the volume of the domain:
print(ad["gas", "volume"].sum())
print(ds.domain_width.prod())

# Now let's make a downsampled plot
plot = yt.SlicePlot(ds, "z", ("gas", "density"), data_source=ad)
plot.save("downsampled.png")
