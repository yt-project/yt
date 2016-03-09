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
import os

from yt.frontends.stream.data_structures import load_particles
from yt.geometry.oct_container import \
    OctreeContainer
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer, \
    ParticleForest
from yt.geometry.oct_container import _ORDER_MAX
from yt.geometry.selection_routines import RegionSelector, AlwaysSelector
from yt.testing import \
    assert_equal, \
    requires_file, \
    assert_true, \
    assert_array_equal
from yt.units.unit_registry import UnitRegistry
from yt.units.yt_array import YTArray
from yt.utilities.lib.geometry_utils import get_morton_indices

import yt.units.dimensions as dimensions
import yt.data_objects.api

NPART = 32**3
DLE = np.array([0.0, 0.0, 0.0])
DRE = np.array([10.0, 10.0, 10.0])
DW = (DRE-DLE)
dx = DW/(2**_ORDER_MAX)

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

    def set_edges(self, file_id, dx = 0.1):
        self.left_edge = YTArray([file_id + dx, 0.0, 0.0],
                                 'code_length', registry=self.ds.unit_registry)
        self.right_edge = YTArray([file_id+1 - dx, self.nfiles, self.nfiles],
                                  'code_length', registry=self.ds.unit_registry)

class FakeBoxRegion:
    def __init__(self, nfiles, DLE, DRE):
        self.ds = FakeDS()
        self.ds.domain_left_edge = YTArray(DLE, "code_length",
                                           registry=self.ds.unit_registry)
        self.ds.domain_right_edge = YTArray(DRE, "code_length",
                                            registry=self.ds.unit_registry)
        self.ds.domain_width = self.ds.domain_right_edge - \
                               self.ds.domain_left_edge
        self.nfiles = nfiles

    def set_edges(self, center, width):
        self.left_edge = self.ds.domain_left_edge + self.ds.domain_width*(center-width/2)
        self.right_edge = self.ds.domain_left_edge + self.ds.domain_width*(center+width/2)

def FakeForest(npart, nfiles, order1, order2, file_order='grid', 
               DLE=None, DRE=None, verbose=False):
    np.random.seed(int(0x4d3d3d3))
    N = (1<<order1)
    if DLE is None: DLE = np.array([0.0, 0.0, 0.0])
    if DRE is None: DRE = np.array([1.0, 1.0, 1.0])
    DW = DRE - DLE
    reg = ParticleForest(DLE, DRE,
                         [N, N, N], nfiles,
                         index_order1 = order1,
                         index_order2 = order2)
    # Create positions for each files
    positions = {}
    if file_order == 'grid':
        # TODO: handle 'remainder' particles
        nYZ = int(np.sqrt(npart/nfiles))
        nR = npart - nYZ*nYZ*nfiles
        div = DW/nYZ
        Y, Z = np.mgrid[DLE[1] + 0.1*div[1] : DRE[1] - 0.1*div[1] : nYZ * 1j,
                        DLE[2] + 0.1*div[2] : DRE[2] - 0.1*div[2] : nYZ * 1j]
        X = 0.5 * div[0] * np.ones(Y.shape, dtype="float64")
        pos = np.array([X.ravel(),Y.ravel(),Z.ravel()],
                       dtype="float64").transpose()
        for i in range(nfiles):
            positions[i] = np.empty_like(pos)
            # positions[i][:,:] = pos
            np.copyto(positions[i],pos)
            pos[:,0] += div[0]
    elif file_order == 'sliced':
        nPF = int(npart/nfiles)
        nR = npart % nfiles
        if verbose: print("{}/{} remainder particles put in first file".format(nR,npart))
        # First file gets extra particles
        pos = np.random.normal(0.5, scale=0.05, size=(nPF+nR,3)) * (DRE-DLE) + DLE
        iLE = DLE[0]
        iRE = iLE + DW[0]/float(nfiles)
        np.clip(pos[:,0], iLE, iRE, pos[:,0])
        for i in range(1,3):
            np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
        positions[0] = pos
        # Other files
        for ifile in range(1,nfiles):
            pos = np.random.normal(0.5, scale=0.05, size=(nPF,3)) * (DRE-DLE) + DLE
            iLE = DLE[0]+DW[0]*float(ifile)/float(nfiles)
            iRE = iLE + DW[0]/float(nfiles)
            np.clip(pos[:,0], iLE, iRE, pos[:,0])
            for i in range(1,3):
                np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
            positions[ifile] = pos
    elif file_order == 'random':
        nPF = int(npart/nfiles)
        nR = npart % nfiles
        if verbose: print("{}/{} remainder particles put in first file".format(nR,npart))
        # First file gets extra particles
        pos = np.random.normal(0.5, scale=0.05, size=(nPF+nR,3)) * (DRE-DLE) + DLE
        for i in range(3):
            np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
        positions[0] = pos
        # Other files
        for ifile in range(nfiles):
            pos = np.random.normal(0.5, scale=0.05, size=(nPF,3)) * (DRE-DLE) + DLE
            for i in range(3):
                np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
            positions[ifile] = pos
    elif file_order == 'octtree': # Not really...
        nfiles_per_dim = int(np.round(nfiles**(1.0/3.0)))
        if (nfiles_per_dim**3) > nfiles:
            nfiles_per_dim -= 1
        nPF = int(npart/(nfiles_per_dim**3))
        nR = npart - nPF*(nfiles_per_dim**3)
        if verbose: print("{}/{} remainder particles put in first file".format(nR,npart))
        div = DW/nfiles_per_dim
        iLE = DLE
        iRE = DLE + div
        ifile = 0
        for xfile in range(nfiles_per_dim):
            iLE[0] = DLE[0] + xfile*div[0]
            iRE[0] = iLE[0] + div[0]
            for yfile in range(nfiles_per_dim):
                iLE[1] = DLE[1] + yfile*div[1]
                iRE[1] = iLE[1] + div[1]
                for zfile in range(nfiles_per_dim):
                    iLE[2] = DLE[2] + zfile*div[2]
                    iRE[2] = iLE[2] + div[2]
                    if ifile == 0:
                        pos = np.random.normal(0.5, scale=0.05, size=(nPF+nR,3)) * (iRE-iLE) + iLE
                    else:
                        pos = np.random.normal(0.5, scale=0.05, size=(nPF,3)) * (iRE-iLE) + iLE
                    for i in range(3):
                        np.clip(pos[:,i], iLE[i], iRE[i], pos[:,i])
                    positions[ifile] = pos
                    ifile += 1
        if verbose: print("Filled {}/{} files".format(ifile,nfiles))
        while ifile < nfiles:
            positions[ifile] = np.zeros((0,3), dtype=pos.dtype)
            ifile += 1
    elif file_order == 'kdtree': # Not really...
        ndense = 1
        ndense_ref = 2
        nfiles_per_dim = int(np.round(nfiles**(1.0/3.0)))
        if (nfiles_per_dim**3) > nfiles:
            nfiles_per_dim -= 1
        nfiles_per_dim -= ndense
        nPF = int(npart/(nfiles_per_dim**3 + (ndense**3)*(ndense_ref**3 - 1)))
        nR = npart - nPF*(nfiles_per_dim**3 + (ndense**3)*(ndense_ref**3 - 1))
        if verbose: print("{}/{} remainder particles put in first file".format(nR,npart))
        div = DW/nfiles_per_dim
        div_dense = div/ndense_ref
        iLE = DLE
        iRE = DLE + div
        iiLE = iLE
        iiRE = iLE + div_dense
        ifile = 0
        for xfile in range(nfiles_per_dim):
            iLE[0] = DLE[0] + xfile*div[0]
            iRE[0] = iLE[0] + div[0]
            for yfile in range(nfiles_per_dim):
                iLE[1] = DLE[1] + yfile*div[1]
                iRE[1] = iLE[1] + div[1]
                for zfile in range(nfiles_per_dim):
                    iLE[2] = DLE[2] + zfile*div[2]
                    iRE[2] = iLE[2] + div[2]
                    if (nfiles_per_dim/2 <= xfile < (nfiles_per_dim/2)+ndense) and \
                       (nfiles_per_dim/2 <= yfile < (nfiles_per_dim/2)+ndense) and \
                       (nfiles_per_dim/2 <= zfile < (nfiles_per_dim/2)+ndense):
                        for xdense in range(ndense_ref):
                            iiLE[0] = iLE[0] + xdense*div_dense[0]
                            iiRE[0] = iiLE[0] + div_dense[0]
                            for ydense in range(ndense_ref):
                                iiLE[1] = iLE[1] + ydense*div_dense[1]
                                iiRE[1] = iiLE[1] + div_dense[1]
                                for zdense in range(ndense_ref):
                                    iiLE[2] = iLE[2] + zdense*div_dense[2]
                                    iiRE[2] = iiLE[2] + div_dense[2]
                                    if ifile == 0:
                                        pos = np.random.normal(0.5, scale=0.05, size=(nPF+nR,3)) * (iiRE-iiLE) + iiLE
                                    else:
                                        pos = np.random.normal(0.5, scale=0.05, size=(nPF,3)) * (iiRE-iiLE) + iiLE
                                    for i in range(3):
                                        np.clip(pos[:,i], iiLE[i], iiRE[i], pos[:,i])
                                    positions[ifile] = pos
                                    ifile += 1
                    else:
                        if ifile == 0:
                            pos = np.random.normal(0.5, scale=0.05, size=(nPF+nR,3)) * (iRE-iLE) + iLE
                        else:
                            pos = np.random.normal(0.5, scale=0.05, size=(nPF,3)) * (iRE-iLE) + iLE
                        for i in range(3):
                            np.clip(pos[:,i], iLE[i], iRE[i], pos[:,i])
                        positions[ifile] = pos
                        ifile += 1
        if verbose: print("Filled {}/{} files".format(ifile,nfiles))
        while ifile < nfiles:
            positions[ifile] = np.zeros((0,3), dtype=pos.dtype)
            ifile += 1
    else:
        raise ValueError("Unsupported value {} for input parameter 'file_order'".format(file_order))
    # Coarse index
    max_npart = positions[0].shape[0]
    sub_mi1 = np.zeros(max_npart, "uint64")
    sub_mi2 = np.zeros(max_npart, "uint64")
    cp = 0
    for i in range(nfiles):
        reg._coarse_index_data_file(positions[i], i)
        cp += positions[i].shape[0]
    if verbose: print("{} particles in total".format(cp))
    reg.find_collisions_coarse(verbose=verbose)
    # Refined index
    for i in range(nfiles):
        reg._refined_index_data_file(positions[i], 
                                     reg.masks.sum(axis=1).astype('uint8'),
                                     sub_mi1, sub_mi2, i)
    reg.find_collisions_refined(verbose=verbose)
    return reg

def time_selection_order2():
    verbose = True
    nfiles = 256
    npart = 128**3
    order1 = 6
    list_order2 = [0,1,2,3,4,5,6,7,8]
    DLE = [0.0, 0.0, 0.0]
    DRE = [1.0, 1.0, 1.0]
    ngz = 1
    # 'grid' also available, but should only be used to test selectors for
    # correctness since it is time consuming
    #file_orders = ['kdtree','octtree']
    file_order = 'random'
    fake_regions = []
    for c,r in [(0.5,0.1),(0.3,0.1),(0.5,0.01),(0.5,0.2),(0.5,0.5),(0.5,1.0)]:
        fr = FakeBoxRegion(nfiles, DLE, DRE)
        fr.set_edges(c,r)
        fake_regions.append(fr)
    if verbose: print("Timing differences due to order2")
    out = {}
    for order2 in list_order2:
        if verbose: print("order2 = {}".format(order2))
        reg = FakeForest(npart, nfiles, order1, order2, file_order=file_order,
                         verbose=verbose)
        t1 = time.time()
        for fr in fake_regions:
            selector = RegionSelector(fr)
            df, gf = reg.identify_data_files(selector, ngz=ngz)
        t2 = time.time()
        out[file_order] = dict(t = t2-t1,
                               ndf = len(df),
                               ngf = len(gf))
        if verbose: print("order2 = {}: {} s, {} files, {} ghost files".format(file_order,t2-t1,len(df),len(gf)))

def time_selection_fileorder():
    verbose = False
    nfiles = 256
    npart = 128**3
    order1 = 6
    order2 = 4
    DLE = [0.0, 0.0, 0.0]
    DRE = [1.0, 1.0, 1.0]
    ngz = 1
    # 'grid' also available, but should only be used to test selectors for
    # correctness since it is time consuming
    #file_orders = ['kdtree','octtree']
    file_orders = ['sliced','random','octtree','kdtree']
    fake_regions = []
    for c,r in [(0.5,0.1),(0.3,0.1),(0.5,0.01),(0.5,0.2),(0.5,0.5),(0.5,1.0)]:
        fr = FakeBoxRegion(nfiles, DLE, DRE)
        fr.set_edges(c,r)
        fake_regions.append(fr)
    if verbose: print("Timing differences due to order of particles.")
    out = {}
    for file_order in file_orders:
        if verbose: print(file_order)
        reg = FakeForest(npart, nfiles, order1, order2, file_order=file_order,
                         verbose=verbose)
        t1 = time.time()
        for fr in fake_regions:
            selector = RegionSelector(fr)
            df, gf = reg.identify_data_files(selector, ngz=ngz)
        t2 = time.time()
        out[file_order] = dict(t = t2-t1,
                               ndf = len(df),
                               ngf = len(gf))
        if verbose: print("{}: {} s, {} files, {} ghost files".format(file_order,t2-t1,len(df),len(gf)))

def test_particle_regions():
    np.random.seed(int(0x4d3d3d3))
    dx = 0.1
    verbose = False
    # We are going to test having 31, 127, 128 and 257 data files.
    # for nfiles in [2, 31, 32, 33, 127, 128, 129]:
    #for nfiles in [2, 31, 32, 33]:
    for nfiles in [2, 31, 127, 128, 129]:
        if verbose: print("nfiles = {}".format(nfiles))
        # Now we create particles 
        # Note: we set order1 to log2(nfiles) here for testing purposes. 
        # Inside the code we set it to min(log2(nfiles), 8)?
        # langmm: this is not strictly true anymore
        # TODO: remove the dims parameter (no longer used by Forest)
        N = nfiles
        order1 = int(np.ceil(np.log2(N))) # Ensures zero collisions
        order2 = 1 # No overlap for N = nfiles
        exact_division = (N == (1 << order1))
        div = float(nfiles)/float(1 << order1)
        reg = FakeForest(nfiles**3, nfiles, order1, order2, file_order='grid',
                         DLE=np.array([0.0, 0.0, 0.0]),
                         DRE=np.array([nfiles, nfiles, nfiles]), 
                         verbose=verbose)
        # Loop over regions selecting single files
        fr = FakeRegion(nfiles)
        for i in range(nfiles):
            fr.set_edges(i, dx)
            selector = RegionSelector(fr)
            df, gf = reg.identify_data_files(selector, ngz=1)
            if exact_division:
                yield assert_equal, len(df), 1, "selector {}, number of files".format(i)
                yield assert_equal, df[0], i, "selector {}, file selected".format(i)
                if i == 0:
                    yield assert_equal, len(gf), 1, "selector {}, number of ghost files".format(i)
                    yield assert_equal, gf[0], i+1, "selector {}, ghost files".format(i)
                elif i == (nfiles - 1):
                    yield assert_equal, len(gf), 1, "selector {}, number of ghost files".format(i)
                    yield assert_equal, gf[0], i-1, "selector {}, ghost files".format(i)
                else:
                    yield assert_equal, len(gf), 2, "selector {}, number of ghost files".format(i)
                    yield assert_equal, gf[0], i-1, "selector {}, ghost files".format(i)
                    yield assert_equal, gf[1], i+1, "selector {}, ghost files".format(i)
            else:
                lf_frac = np.floor(float(fr.left_edge[0])/div)*div
                rf_frac = np.floor(float(fr.right_edge[0])/div)*div
                # Selected files
                lf = int(np.floor(lf_frac) if ((lf_frac % 0.5) == 0) else np.round(lf_frac))
                rf = int(np.floor(rf_frac) if ((rf_frac % 0.5) == 0) else np.round(rf_frac))
                if (rf+0.5) >= (rf_frac+div): rf -= 1
                if (lf+0.5) <= (lf_frac-div): lf += 1
                df_ans = np.arange(max(lf,0),min(rf+1,nfiles))
                # print df, df_ans
                # print lf_frac, lf, rf_frac, rf, lf_frac-div, (rf_frac+div)
                yield assert_array_equal, df, df_ans, "selector {}, file array".format(i)
                # Ghost zones selected files
                lf_ghost = int(max(np.floor(lf_frac - div) if (((lf_frac-div) % 0.5) == 0) else np.round(lf_frac - div),0))
                rf_ghost = int(min(np.floor(rf_frac + div) if (((rf_frac+div) % 0.5) == 0) else np.round(rf_frac + div),nfiles-1))
                if (rf_ghost+0.5) >= (rf_frac+2*div): rf_ghost -= 1
                gf_ans = []
                if lf_ghost < lf: gf_ans.append(lf_ghost)
                if rf_ghost > rf: gf_ans.append(rf_ghost)
                gf_ans = np.array(gf_ans)
                yield assert_array_equal, gf, gf_ans, "selector {}, ghost file array".format(i)

        # print reg.masks.shape
        # for mask in reg.masks:
        #     print mask.shape
        #     maxs = np.unique(mask.max(axis=-1).max(axis=-1))
        #     mins = np.unique(mask.min(axis=-1).min(axis=-1))
        #     yield assert_equal, maxs, mins
        #     yield assert_equal, maxs, np.unique(mask)

def test_save_load_bitmap():
    verbose = False
    fname_fmt = "temp_bitmasks{}.dat"
    i = 0
    fname = fname_fmt.format(i)
    while os.path.isfile(fname):
        i += 1
        fname = fname_fmt.format(i)
    np.random.seed(int(0x4d3d3d3))
    nfiles = 32
    order1 = 2
    order2 = 2 # Maximum collisions
    pos = np.random.normal(0.5, scale=0.05, size=(NPART/nfiles,3)) * (DRE-DLE) + DLE
    pos[:,0] = (DW[0]/nfiles)/2
    for i in range(3):
        np.clip(pos[:,i], DLE[i], DRE[i], pos[:,i])
    reg0 = ParticleForest(DLE, DRE, 3*[1<<order1], nfiles,
                          index_order1 = order1,
                          index_order2 = order2)
    # Coarse index
    for i in range(nfiles):
        reg0._coarse_index_data_file(pos, i)
        pos[:,0] += (DW[0]/nfiles)
    pos[:,0] = (DW[0]/nfiles)/2
    reg0.find_collisions_coarse(verbose=verbose)
    # Refined index
    max_npart = pos.shape[0]
    sub_mi1 = np.zeros(max_npart, "uint64")
    sub_mi2 = np.zeros(max_npart, "uint64")
    for i in range(0,nfiles):
        reg0._refined_index_data_file(pos, reg0.masks.sum(axis=1).astype('uint8'),
                                      sub_mi1, sub_mi2, i)
        pos[:,0] += (DW[0]/nfiles)
    pos[:,0] = (DW[0]/nfiles)/2
    reg0.find_collisions_refined(verbose=verbose)
    # Save
    reg0.save_bitmasks(fname)
    # Load
    reg1 = ParticleForest(DLE, DRE, 3*[1<<order1], nfiles,
                          index_order1 = order1,
                          index_order2 = order2)
    reg1.load_bitmasks(fname)
    # Check equality
    yield assert_true, reg0.iseq_bitmask(reg1)
    # Remove file
    os.remove(fname)

if __name__=="__main__":
    for i in test_add_particles_random():
        i[0](*i[1:])
    time.sleep(1)

os33 = "snapshot_033/snap_033.0.hdf5"
@requires_file(os33)
def test_get_smallest_dx():
    ds = yt.load(os33)
    yield assert_equal, ds.index.get_smallest_dx(), \
        ds.domain_width / (ds.domain_dimensions*2.**(ds.index.max_level))

# TODO: Change these!!!!
bc94 = "/mnt/gv0/mturk/big_cosmo/snapdir_094/snap_lcdma_1024_094.0"
bc94_coll = "/root/projects/bitmap/big_cosmo_bitmask_7_5_coll.dat"

@requires_file(bc94)
@requires_file(bc94_coll)
def test_initialize_index():
    order1 = 7
    order2 = 5
    ds = yt.GadgetDataset(bc94, long_ids = True)
    ds.index._initialize_index(order1=order1, order2=order2)
    reg1 = ds.index.regions
    N = ds.index.ds.domain_dimensions / (1<<ds.index.ds.over_refine_factor)
    reg0 = ParticleForest(ds.domain_left_edge, ds.domain_right_edge,
                          N, len(ds.index.data_files), ds.over_refine_factor,
                          ds.n_ref, index_order1=order1, index_order2=order2)
    reg0.load_bitmasks(fname=bc94_coll)
    yield assert_true, reg0.iseq_bitmask(reg1)

# To avoid loading
class FakeBC94DS:
    unit_registry = UnitRegistry()
    unit_registry.add('code_length', 1.0, dimensions.length)
    domain_left_edge = YTArray([0.0, 0.0, 0.0], "code_length",
                               registry=unit_registry)
    domain_right_edge = YTArray([135.54, 135.54, 135.54], "code_length",
                                registry=unit_registry)
    domain_width = YTArray([135.54, 135.54, 135.54], "code_length",
                           registry=unit_registry)
    domain_center = YTArray([ 67.77,  67.77,  67.77], "code_length",
                            registry=unit_registry)
    periodicity = (True, True, True)
    over_refine_factor = 1
    n_ref = 64
    nfiles = 512
    order1 = 7
    order2 = 5
    domain_dimensions = np.array(3*[nfiles], dtype='int32')
    default_fluid_type = 'gas'

# class FakeBC94SphericalRegion:
#     #from yt.geometry.selection_routines import SphereSelector
#     from yt.geometry.selection_routines import sphere_selector
#     def __init__(self, c, r):
#         self.ds = FakeBC94DS()
#         self.nfiles = self.ds.nfiles
#         self.center = self.ds.domain_center + c*self.ds.domain_width
#         self.radius = r*self.ds.domain_width[0]
#         #self.selector = SphereSelector(self)
#         self.selector = sphere_selector(self)

# TODO: remove dependence on bc94 (only use bitmask)
@requires_file(bc94)
@requires_file(bc94_coll)
def test_fill_masks():
    from yt.utilities.lib.ewah_bool_wrap import BoolArrayCollection
    from yt.geometry.particle_oct_container import ParticleForestSelector
    from yt.data_objects.selection_data_containers import YTSphere
    ds_empty = FakeBC94DS()

    order1 = 7
    order2 = 5
    ngz = 1
    ds = yt.GadgetDataset(bc94, long_ids = True)
    ds.index._initialize_index(fname=bc94_coll, order1=order1, order2=order2)
    print "default_fluid_type",getattr(ds,"default_fluid_type")
    N = ds_empty.domain_dimensions/(1<<ds_empty.over_refine_factor)
    reg = ds.index.regions
    reg0 = ParticleForest(ds_empty.domain_left_edge, ds_empty.domain_right_edge,
                          N, ds_empty.nfiles, ds_empty.over_refine_factor,
                          ds_empty.n_ref, index_order1=ds_empty.order1, index_order2=ds_empty.order2)
    reg0.load_bitmasks(fname=bc94_coll)
    tests_sph = {(0,0.9/float(1 << order1)):  8,
                 (0,1.0/float(1 << order1)): 20,
                 (0,1.1/float(1 << order1)): 32,
                 (0.5/float(1 << order1),0.49/float(1 << order1)): 1,
                 (0.5/float(1 << order1),0.50/float(1 << order1)): 4,
                 (0.5/float(1 << order1),0.51/float(1 << order1)): 7,
                 (0,1.9/float(1 << order1)): 64,
                 # (0,2.0/float(1 << order1)): 76, # floating point equality...
                 (0,2.1/float(1 << order1)): 88,
                 (0.5/float(1 << order1),1.49/float(1 << order1)): 27,
                 (0.5/float(1 << order1),1.50/float(1 << order1)): 30,
                 (0.5/float(1 << order1),1.51/float(1 << order1)): 33}
    for (c, r), nc_s in tests_sph.items():
        mm_s = BoolArrayCollection()
        mm_g = BoolArrayCollection()
        mm_s0 = BoolArrayCollection()
        mm_g0 = BoolArrayCollection()
        center = ds.domain_center + c*ds.domain_width
        radius = r*ds.domain_width[0]
        sp = ds.sphere(center, radius)
        #sp0 = FakeBC94SphericalRegion(c, r)
        sp0 = YTSphere(center, radius, ds=ds_empty)
        ms = ParticleForestSelector(sp.selector, reg, ngz=ngz)
        ms0 = ParticleForestSelector(sp0.selector, reg0, ngz=ngz)
        ms.fill_masks(mm_s, mm_g)
        ms0.fill_masks(mm_s0, mm_g0)
        yield assert_equal, mm_s.count_coarse(), nc_s
        print(c,r,nc_s,mm_s.count_coarse(),"succeeded")
        yield assert_equal, mm_s0.count_coarse(), nc_s
        print(c,r,nc_s,mm_s0.count_coarse(),"succeeded")
