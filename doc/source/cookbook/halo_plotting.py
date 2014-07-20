import yt
from yt.analysis_modules.halo_analysis.halo_catalog import HaloCatalog

# Load the dataset
ds = yt.load("Enzo_64/RD0006/RedshiftOutput0006")

# Load the halo list from a rockstar output for this dataset
halos = yt.load('rockstar_halos/halos_0.0.bin')

# Create the halo catalog from this halo list
hc = HaloCatalog(halos_ds=halos)
hc.load()

# Create a projection with the halos overplot on top
p = yt.ProjectionPlot(ds, "x", "density")
p.annotate_halos(hc)
p.save()
