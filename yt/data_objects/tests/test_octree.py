"""
Tests for YTOctree


"""

# ----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from yt.testing import \
    fake_sph_grid_ds, fake_sph_orientation_ds, assert_equal
from yt.data_objects.construction_data_containers import YTOctree
import tempfile
import os
import shutil
import numpy as np

def test_building_tree():
    '''
    Test function to build an octree and make sure correct number of particles
    '''
    ds = fake_sph_grid_ds()
    print(ds.field_list)
    octree = ds.octree(n_ref=1)
    assert(type(octree) == YTOctree)
    assert(octree['x'].shape[0] == 456)

def test_saving_loading():
    '''
    This builds an octree, writes to file, reloads and ensure that the reloaded
    octree is the same as the initial built tree.
    '''
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = fake_sph_grid_ds()
    ds.tree_filename = tmpdir+"test.octree"
    ds._file_hash = 1
    octree = ds.octree(n_ref=1)

    ds2 = fake_sph_grid_ds()
    ds2.tree_filename = tmpdir+"test.octree"
    ds2._file_hash = 1
    octree_loaded = ds2.octree(n_ref=1)

    assert(octree == octree_loaded)
    assert(octree_loaded.loaded)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

def test_sph_interpolation_scatter():
    '''
    Just generate an octree, perform some SPH interpolation and check with some
    answer testing
    '''

    ds = fake_sph_grid_ds(hsml_factor=15.0)
    ds.use_sph_normalization = False
    octree = ds.octree(n_ref=1, over_refine_factor=1)
    density = octree[('all', 'density')]
    print(density)

def test_sph_interpolation_gather():
    '''
    Just generate an octree, perform some SPH interpolation and check with some
    answer testing
    '''
    ds = fake_sph_grid_ds(hsml_factor=15.0)
    ds.sph_smoothing_style = 'gather'
    ds.num_neighbors = 5
    ds.use_sph_normalization = False
    octree = ds.octree(n_ref=1, over_refine_factor=1)
    density = octree[('all', 'density')]
    print(density)

def test_over_refine_factor():
    '''
    Ensure that the octree over refine factor is behaving as expected
    '''
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=1, over_refine_factor=2)
    num_cells = octree['x'].shape[0]
    assert(num_cells == 3648)

def test_density_factor():
    '''
    Ensure the dense tree functionality is working
    '''
    ds = fake_sph_grid_ds()
    octree = ds.octree(n_ref=1, density_factor=2)
    num_cells = octree['x'].shape[0]
    print(num_cells)
    assert(num_cells == 1)
