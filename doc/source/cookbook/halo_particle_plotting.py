"""
This is a simple mechanism for overplotting the particles belonging only to
halos.  For more information, see :ref:`halo_finding`.
"""
from yt.mods import * # set up our namespace

pf = load("Enzo_64/DD0043/data0043")

halos = HaloFinder(pf)

p = ProjectionPlot(pf, "x", "density")
p.annotate_hop_circles(halos)
p.annotate_hop_particles(halos, max_number=100)
p.save()
