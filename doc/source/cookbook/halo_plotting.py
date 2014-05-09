"""
This is a mechanism for plotting circles representing identified particle halos
on an image.  For more information, see :ref:`halo_finding`.
"""
from yt.mods import * # set up our namespace

data_pf = load("Enzo_64/RD0006/RedshiftOutput0006")

halo_pf = load('rockstar_halos/halos_0.0.bin')

hc - HaloCatalog(halos_pf = halo_pf)
hc.load()

p = ProjectionPlot(pf, "x", "density")
p.annotate_halos(hc)
p.save()
