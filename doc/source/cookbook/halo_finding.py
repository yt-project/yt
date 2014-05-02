"""
This script shows the simplest way of getting halo information.  For more
information, see :ref:`halo_finding`.
"""
import yt

ds = yt.load("Enzo_64/DD0043/data0043")

halos = yt.HaloFinder(ds)
halos.write_out("%s_halos.txt" % ds)
