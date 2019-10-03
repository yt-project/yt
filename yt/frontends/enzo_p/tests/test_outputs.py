"""
Title: test_enzo_p.py
Purpose: Enzo-P frontend tests
    Copyright (c) 2017, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import numpy as np
import pytest

from yt.utilities.on_demand_imports import \
    _h5py as h5py

from yt.testing import \
    assert_equal, \
    requires_file, \
    assert_array_equal
from yt.frontends.enzo_p.api import EnzoPDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
hello_world = "hello-0210/hello-0210.block_list"
ep_cosmo = "ENZOP_DD0140/ENZOP_DD0140.block_list"


# Global field info
_fields = ("density", "total_energy",
           "velocity_x", "velocity_y")
_pfields = ("particle_position_x", "particle_position_y",
            "particle_position_z", "particle_velocity_x",
            "particle_velocity_y", "particle_velocity_z")


#============================================
#                TestEnzoP
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestEnzoP(fw.AnswerTest):
    #-----
    # test_EnzoDataset
    #-----
    @requires_file(hello_world)
    def test_EnzoPDataset(self, ds_hello_world):
        assert isinstance(ds_hello_world, EnzoPDataset)

    #-----
    # test_hello_world
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hello_world)
    def test_hello_world(self, ds_hello_world):
        self.hashes['pixelized_projection_values'] = OrderedDict()
        self.hashes['field_values'] = OrderedDict()
        ds = ds_hello_world
        dso = [ None, ("sphere", ("max", (0.25, 'unitary')))]
        for d in dso:
            self.hashes['pixelized_projection_values'][d] = OrderedDict()
            self.hashes['field_values'][d] = OrderedDict()
            for f in _fields:
                self.hashes['pixelized_projection_values'][d][f] = OrderedDict()
                for a in [0, 1, 2]:
                    self.hashes['pixelized_projection_values'][d][f][a] = OrderedDict()
                    for w in [None, "density"]:
                        ppv_hd = self.pixelized_projection_values_test(ds, a, f, w, d)
                    self.hashes['pixelized_projection_values'][d][f][a][w] = ppv_hd
                fv_hd = self.field_values_test(ds, f, d)
                self.hashes['field_values'][d][f] = fv_hd
            dobj = utils.create_obj(ds, d)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert_equal(s1, s2)

    #-----
    # test_particle_fields
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(ep_cosmo)
    def test_particle_fields(self, ds_ep_cosmo):
        self.hashes['field_values'] = OrderedDict()
        ds = ds_ep_cosmo
        dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
        for dobj_name in dso:
            self.hashes['field_values'][dobj_name] = OrderedDict()
            for field in _pfields:
                fv_hd = self.field_values_test(ds, field, dobj_name, particle_type=True)
                self.hashes['field_values'][dobj_name][field] = fv_hd
            dobj = utils.create_obj(ds, dobj_name)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert_equal(s1, s2)

    #-----
    # test_hierarchy
    #-----
    @requires_file(hello_world)
    def test_hierarchy(self, ds_hello_world):
        ds = ds_hello_world
        fh = h5py.File(ds.index.grids[0].filename, "r")
        for grid in ds.index.grids:
            assert_array_equal(
                grid.LeftEdge.d, fh[grid.block_name].attrs["enzo_GridLeftEdge"])
            assert_array_equal(
                ds.index.grid_left_edge[grid.id], grid.LeftEdge)
            assert_array_equal(
                ds.index.grid_right_edge[grid.id], grid.RightEdge)
            for child in grid.Children:
                assert (child.LeftEdge >= grid.LeftEdge).all()
                assert (child.RightEdge <= grid.RightEdge).all()
                assert_equal(child.Parent.id, grid.id)
        fh.close()

    #-----
    # test_critical_density
    #-----
    @requires_file(ep_cosmo)
    def test_critical_density(self, ds_ep_cosmo):
        ds = ds_ep_cosmo
        c1 = (ds.r["dark", "particle_mass"].sum() +
              ds.r["gas", "cell_mass"].sum()) / \
              ds.domain_width.prod() / ds.critical_density
        c2 = ds.omega_matter * (1 + ds.current_redshift)**3 / \
          (ds.omega_matter * (1 + ds.current_redshift)**3 + ds.omega_lambda)
        assert np.abs(c1 - c2) / max(c1, c2) < 1e-3
