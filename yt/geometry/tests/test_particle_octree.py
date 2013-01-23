from yt.testing import *
import numpy as np
from yt.geometry.oct_container import ParticleOctreeContainer
import time, os

NPART = 32**3
NDIM = 64
DLE = np.array([0.0, 0.0, 0.0])
DRE = np.array([10.0, 10.0, 10.0])

def test_add_particles_random():
    np.random.seed(int(0x4d3d3d3))
    pos = np.random.normal(0.5, scale=0.05, size=(NPART,3)) * (DRE-DLE) + DLE
    for i in range(3):
        np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
    for ndom in [1, 2, 4, 8]:
        octree = ParticleOctreeContainer((NDIM, NDIM, NDIM), DLE, DRE)
        for dom in range(ndom):
            octree.add(pos[dom::ndom,:], dom)
        octree.finalize()
        # This visits every oct.
        lin_count = octree.linearly_count()
        tc = octree.recursively_count()
        total_count = np.zeros(len(tc), dtype="int32")
        for i in sorted(tc):
            total_count[i] = tc[i]
        yield assert_equal, lin_count, total_count.sum()
        mask = np.ones((total_count.sum(), 8), dtype="bool")
        # This visits every cell -- including those covered by octs.
        level_count  = octree.count_levels(total_count.size-1, -1, mask)
        for dom in range(ndom):
            level_count += octree.count_levels(total_count.size-1, dom, mask)
        yield assert_equal, level_count[0], NDIM**3 * 8
        yield assert_equal, level_count, total_count * 8

if __name__=="__main__":
    for i in test_add_particles_random():
        i[0](*i[1:])
    time.sleep(1)
