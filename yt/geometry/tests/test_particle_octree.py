"""
Tests for particle octree



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
import time

from yt.frontends.stream.data_structures import load_particles
from yt.geometry.oct_container import \
    OctreeContainer
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer, \
    ParticleRegions
from yt.geometry.oct_container import _ORDER_MAX
from yt.geometry.selection_routines import RegionSelector, AlwaysSelector
from yt.testing import \
    assert_equal, \
    requires_file
from yt.units.unit_registry import UnitRegistry
from yt.units.yt_array import YTArray
from yt.utilities.lib.geometry_utils import get_morton_indices

import yt.units.dimensions as dimensions
import yt.data_objects.api

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

def test_save_load_octree():
    np.random.seed(int(0x4d3d3d3))
    pos = np.random.normal(0.5, scale=0.05, size=(NPART,3)) * (DRE-DLE) + DLE
    octree = ParticleOctreeContainer((1, 1, 1), DLE, DRE)
    octree.n_ref = 32
    for i in range(3):
        np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
    # Convert to integers
    pos = np.floor((pos - DLE)/dx).astype("uint64")
    morton = get_morton_indices(pos)
    morton.sort()
    octree.add(morton)
    octree.finalize()
    saved = octree.save_octree()
    loaded = OctreeContainer.load_octree(saved)
    always = AlwaysSelector(None)
    ir1 = octree.ires(always)
    ir2 = loaded.ires(always)
    yield assert_equal, ir1, ir2

    fc1 = octree.fcoords(always)
    fc2 = loaded.fcoords(always)
    yield assert_equal, fc1, fc2

    fw1 = octree.fwidth(always)
    fw2 = loaded.fwidth(always)
    yield assert_equal, fw1, fw2

def test_particle_octree_counts():
    np.random.seed(int(0x4d3d3d3))
    # Eight times as many!
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
        ds = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref)
        dd = ds.all_data()
        bi = dd["io","mesh_id"]
        v = np.bincount(bi.astype("intp"))
        yield assert_equal, v.max() <= n_ref, True
        bi2 = dd["all","mesh_id"]
        yield assert_equal, bi, bi2

def test_particle_overrefine():
    np.random.seed(int(0x4d3d3d3))
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
        ds1 = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref)
        dd1 = ds1.all_data()
        v1 = dict((a, getattr(dd1, a)) for a in _attrs)
        cv1 = dd1["cell_volume"].sum(dtype="float64")
        for over_refine in [1, 2, 3]:
            f = 1 << (3*(over_refine-1))
            ds2 = load_particles(data, 1.0, bbox = bbox, n_ref = n_ref,
                                over_refine_factor = over_refine)
            dd2 = ds2.all_data()
            v2 = dict((a, getattr(dd2, a)) for a in _attrs)
            for a in sorted(v1):
                yield assert_equal, v1[a].size * f, v2[a].size
            cv2 = dd2["cell_volume"].sum(dtype="float64")
            yield assert_equal, cv1, cv2

index_ptype_snap = "snapshot_033/snap_033.0.hdf5"
@requires_file(index_ptype_snap)
def test_particle_index_ptype():
    ds = yt.load(index_ptype_snap)
    ds_all = yt.load(index_ptype_snap, index_ptype="all")
    ds_pt0 = yt.load(index_ptype_snap, index_ptype="PartType0")
    dd = ds.all_data()
    dd_all = ds_all.all_data()
    dd_pt0 = ds_pt0.all_data()
    cv = dd["cell_volume"]
    cv_all = dd_all["cell_volume"]
    cv_pt0 = dd_pt0["cell_volume"]
    yield assert_equal, cv.shape, cv_all.shape
    yield assert_equal, cv.sum(dtype="float64"), cv_pt0.sum(dtype="float64")

class FakeDS:
    domain_left_edge = None
    domain_right_edge = None
    domain_width = None
    unit_registry = UnitRegistry()
    unit_registry.add('code_length', 1.0, dimensions.length)
    periodicity = (False, False, False)

class FakeRegion:
    def __init__(self, nfiles):
        self.ds = FakeDS()
        self.ds.domain_left_edge = YTArray([0.0, 0.0, 0.0], "code_length",
                                           registry=self.ds.unit_registry)
        self.ds.domain_right_edge = YTArray([nfiles, nfiles, nfiles], "code_length",
                                            registry=self.ds.unit_registry)
        self.ds.domain_width = self.ds.domain_right_edge - \
                               self.ds.domain_left_edge
        self.nfiles = nfiles

    def set_edges(self, file_id):
        self.left_edge = YTArray([file_id + 0.1, 0.0, 0.0],
                                 'code_length', registry=self.ds.unit_registry)
        self.right_edge = YTArray([file_id+1 - 0.1, self.nfiles, self.nfiles],
                                  'code_length', registry=self.ds.unit_registry)

def test_particle_regions():
    np.random.seed(int(0x4d3d3d3))
    # We are going to test having 31, 127, 128 and 257 data files.
    for nfiles in [2, 31, 127, 128, 129]:
        # Now we create particles 
        # Note: we set N to nfiles here for testing purposes.  Inside the code 
        # we set it to min(N, 256)
        N = nfiles
        reg = ParticleRegions([0.0, 0.0, 0.0, 0.0],
                              [nfiles, nfiles, nfiles],
                              [N, N, N], nfiles)
        Y, Z = np.mgrid[0.1 : nfiles - 0.1 : nfiles * 1j,
                        0.1 : nfiles - 0.1 : nfiles * 1j]
        X = 0.5 * np.ones(Y.shape, dtype="float64")
        pos = np.array([X.ravel(),Y.ravel(),Z.ravel()],
            dtype="float64").transpose()
        for i in range(nfiles):
            reg.add_data_file(pos, i)
            pos[:,0] += 1.0
        pos[:,0] = 0.5
        fr = FakeRegion(nfiles)
        for i in range(nfiles):
            fr.set_edges(i)
            selector = RegionSelector(fr)
            df = reg.identify_data_files(selector)
            yield assert_equal, len(df), 1
            yield assert_equal, df[0], i
            pos[:,0] += 1.0

        for mask in reg.masks:
            maxs = np.unique(mask.max(axis=-1).max(axis=-1))
            mins = np.unique(mask.min(axis=-1).min(axis=-1))
            yield assert_equal, maxs, mins
            yield assert_equal, maxs, np.unique(mask)

if __name__=="__main__":
    for i in test_add_particles_random():
        i[0](*i[1:])
    time.sleep(1)

def test_position_location():
    np.random.seed(int(0x4d3d3d3))
    pos = np.random.normal(0.5, scale=0.05, size=(NPART,3)) * (DRE-DLE) + DLE
    # Now convert to integers
    data = {}
    bbox = []
    for i, ax in enumerate('xyz'):
        np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
        bbox.append([DLE[i], DRE[i]])
        data["particle_position_%s" % ax] = pos[:,i]
    bbox = np.array(bbox)
    ds = load_particles(data, 1.0, bbox = bbox, over_refine_factor = 2)
    oct_id, all_octs = ds.index.oct_handler.locate_positions(pos)
    for oi in sorted(all_octs):
        this_oct = pos[oct_id == oi]
        assert(np.all(this_oct >= all_octs[oi]["left_edge"]))
        assert(np.all(this_oct <= all_octs[oi]["right_edge"]))

os33 = "snapshot_033/snap_033.0.hdf5"
@requires_file(os33)
def test_get_smallest_dx():
    ds = yt.load(os33)
    yield assert_equal, ds.index.get_smallest_dx(), \
        ds.domain_width / (ds.domain_dimensions*2.**(ds.index.max_level))
