from yt.testing import *
import numpy as np
from yt.geometry.particle_oct_container import ParticleOctreeContainer
from yt.geometry.oct_container import _ORDER_MAX
from yt.utilities.lib.geometry_utils import get_morton_indices
from yt.frontends.stream.api import load_particles
import yt.data_objects.api
import time, os

NPART = 32**3
DLE = np.array([0.0, 0.0, 0.0])
DRE = np.array([10.0, 10.0, 10.0])
dx = (DRE-DLE)/(2**_ORDER_MAX)

def test_add_particles_random():
    np.random.seed(int(0x4d3d3d3))
    pos = np.random.normal(0.5, scale=0.05, size=(NPART,3)) * (DRE-DLE) + DLE
    # Now convert to integers
    for i in range(3):
        np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
    # Convert to integers
    pos = np.floor((pos - DLE)/dx).astype("uint64")
    morton = get_morton_indices(pos)
    morton.sort()
    for ndom in [1, 2, 4, 8]:
        octree = ParticleOctreeContainer((1, 1, 1), DLE, DRE)
        octree.n_ref = 32
        for dom, split in enumerate(np.array_split(morton, ndom)):
            octree.add(split)
        octree.finalize()
        # This visits every oct.
        tc = octree.recursively_count()
        total_count = np.zeros(len(tc), dtype="int32")
        for i in sorted(tc):
            total_count[i] = tc[i]
        yield assert_equal, octree.nocts, total_count.sum()
        # This visits every cell -- including those covered by octs.
        #for dom in range(ndom):
        #    level_count += octree.count_levels(total_count.size-1, dom, mask)
        yield assert_equal, total_count, [1, 8, 64, 64, 256, 536, 1856, 1672]

def test_particle_octree_counts():
    np.random.seed(int(0x4d3d3d3))
    # Eight times as many!
    pos = []
    data = {}
    bbox = []
    for i, ax in enumerate('xyz'):
        DW = DRE[i] - DLE[i]
        LE = DLE[i]
        data["particle_position_%s" % ax] = \
            np.random.normal(0.5, scale=0.05, size=(NPART*8)) * DW + LE
        bbox.append( [DLE[i], DRE[i]] )
    bbox = np.array(bbox)
    for n_ref in [16, 32, 64, 512, 1024]:
        pf = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref)
        dd = pf.h.all_data()
        bi = dd["all","mesh_id"]
        v = np.bincount(bi.astype("int64"))
        yield assert_equal, v.max() <= n_ref, True

def test_particle_overrefine():
    np.random.seed(int(0x4d3d3d3))
    pos = []
    data = {}
    bbox = []
    for i, ax in enumerate('xyz'):
        DW = DRE[i] - DLE[i]
        LE = DLE[i]
        data["particle_position_%s" % ax] = \
            np.random.normal(0.5, scale=0.05, size=(NPART)) * DW + LE
        bbox.append( [DLE[i], DRE[i]] )
    bbox = np.array(bbox)
    _attrs = ('icoords', 'fcoords', 'fwidth', 'ires')
    for n_ref in [16, 32, 64, 512, 1024]:
        pf1 = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref)
        dd1 = pf1.h.all_data()
        v1 = dict((a, getattr(dd1, a)) for a in _attrs)
        cv1 = dd1["CellVolumeCode"].sum(dtype="float64")
        for over_refine in [1, 2, 3]:
            f = 1 << (3*(over_refine-1))
            pf2 = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref,
                                over_refine_factor = over_refine)
            dd2 = pf2.h.all_data()
            v2 = dict((a, getattr(dd2, a)) for a in _attrs)
            for a in sorted(v1):
                yield assert_equal, v1[a].size * f, v2[a].size
            cv2 = dd2["CellVolumeCode"].sum(dtype="float64")
            yield assert_equal, cv1, cv2

if __name__=="__main__":
    for i in test_add_particles_random():
        i[0](*i[1:])
    time.sleep(1)
