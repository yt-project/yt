"""
Enzo-P frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import \
    _h5py as h5py

from yt.testing import \
    assert_equal, \
    requires_file, \
    assert_array_equal
from yt.utilities.answer_testing.framework import \
    create_obj, \
    data_dir_load, \
    requires_ds, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest
from yt.frontends.enzo_p.api import EnzoPDataset

_fields = ("density", "total_energy",
           "velocity_x", "velocity_y")

hello_world = "hello-0200/hello-0200.block_list"


@requires_file(hello_world)
def test_EnzoPDataset():
    assert isinstance(data_dir_load(hello_world), EnzoPDataset)

@requires_ds(hello_world)
def test_hello_world():
    ds = data_dir_load(hello_world)

    dso = [ None, ("sphere", ("max", (0.25, 'unitary')))]
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        hello_world, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(hello_world, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)

@requires_file(hello_world)
def test_hierarchy():
    ds = data_dir_load(hello_world)

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
