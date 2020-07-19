import yt

# Load the dataset
ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")

# Load the halo list from a rockstar output for this dataset
halos = yt.load("rockstar_halos/halos_0.0.bin")

# Create a projection with the halos overplot on top
p = yt.ProjectionPlot(ds, "x", "density")
p.annotate_halos(halos)
p.save()
