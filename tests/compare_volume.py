## This script can be deleted once grid_visitors is ready to be merged in.
## It is set up to compare the way the chunking systems work.
##
## There are two things we want to check:
##   1) Do we get the same values (potentially in different order) regardless
##      of how we chunk the data?
##   2) Do we get the same values in the same order regardless of how we chunk
##      the data?
## We do this by comparing three specific things:
##   - The sum of mass; this requires volume and density be in the same order.
##   - The unsorted morton indices (which should be unique)
##   - The sorted morton indices
## There are three chunking systems we need to study:
##   - io
##   - spatial
##   - all

import numpy as np
import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
dd = ds.all_data()

mi1 = dd["morton_index"]
mi2 = []
mi3 = []

m1 = dd["cell_mass"].sum()
m2 = 0.0 * m1
m3 = 0.0 * m1

for c in dd.chunks([], "io"):
    mi2.append(c["morton_index"])
    m2 += c["cell_mass"].sum()

for c in dd.chunks([], "spatial"):
    mi3.append(c["morton_index"][c._current_chunk.objs[0].child_mask])

# Presently disabled, as there are issues with spatial chunking and
# double-iterating.  This is the next problem to be solved.
for i, c in enumerate(dd.chunks([], "spatial")):
    continue
    print("Chunk ...", i)
    m3 += c["cell_mass"][c._current_chunk.objs[0].child_mask].sum()

print("Masses all-io-sp", m1, m2, m3)

mi2 = yt.uconcatenate(mi2)
mi3 = yt.uconcatenate(mi3)

print("Morton all-io", (mi1 == mi2).all())
print("Morton all-sp", (mi1 == mi3).all())
print("Morton io-sp ", (mi2 == mi3).all())
mi1.sort()
mi2.sort()
mi3.sort()
print("Morton sorted all-io", (mi1 == mi2).all())
print("Morton sorted all-sp", (mi1 == mi3).all())
