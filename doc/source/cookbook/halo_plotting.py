"""
This is a mechanism for plotting circles representing identified particle halos
on an image.  For more information, see :ref:`halo_finding`.
"""
from yt.mods import * # set up our namespace

pf = load("Enzo_64/DD0043/data0043")

halos = HaloFinder(pf)

p = ProjectionPlot(pf, "z", "density")
p.annotate_hop_circles(halos)
p.save()
