"""
This script shows the simplest way of getting halo information.  For more
information, see :ref:`halo_finding`.
"""
from yt.mods import * # set up our namespace

pf = load("Enzo_64/DD0043/data0043")

halos = HaloFinder(pf)
halos.write_out("%s_halos.txt" % pf)
