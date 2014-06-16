"""
This is a mechanism for plotting circles representing identified particle halos
on an image.  For more information, see :ref:`halo_finding`.
"""
from yt.mods import * # set up our namespace

data_ds = load("Enzo_64/RD0006/RedshiftOutput0006")

halo_ds = load('rockstar_halos/halos_0.0.bin')

hc - HaloCatalog(halos_ds = halo_ds)
hc.load()

p = ProjectionPlot(ds, "x", "density")
p.annotate_halos(hc)
p.save()
