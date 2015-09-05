"""
Gadget frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import requires_file, assert_equal
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.gadget.api import GadgetHDF5Dataset, GadgetDataset

isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"
isothermal_bin = "IsothermalCollapse/snap_505"
gdg = "GadgetDiskGalaxy/snapshot_0200.hdf5"

_fields = (
    ("gas", "density")
    ("gas", "temperature"),
    ('gas', 'velocity_magnitude'),
    ("deposit", "all_density"),
    ("deposit", "all_count"),
    ("deposit", "all_cic")
    ("deposit", "PartType0_density"),
    ("deposit", "PartType4_density")
)

iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])
gdg_kwargs = dict(bounding_box=[[-1e5, 1e5], [-1e5, 1e5], [-1e5, 1e5]])


@requires_file(isothermal_h5)
@requires_file(isothermal_bin)
def test_GadgetDataset():
    assert isinstance(data_dir_load(isothermal_h5, kwargs=iso_kwargs),
                      GadgetHDF5Dataset)
    assert isinstance(data_dir_load(isothermal_bin, kwargs=iso_kwargs),
                      GadgetDataset)


@requires_ds(isothermal_h5)
def test_iso_collapse():
    ds = data_dir_load(isothermal_h5, kwargs=iso_kwargs)
    yield assert_equal, str(ds), "snap_505"
    dso = [None, ("sphere", ("c", (0.1, 'unitary')))]
    dd = ds.all_data()
    yield assert_equal, dd["particle_position"].shape, (2**17, 3)
    tot = sum(dd[ptype, "particle_position"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    yield assert_equal, tot, (2**17)
    for dobj_name in dso:
        for field in _fields:
            if 'PartType4' in field[1]:
                continue
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        isothermal_h5, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(isothermal_h5, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2


@requires_ds(gdg, big_data=True)
def test_snapshot_200():
    ds = data_dir_load(gdg)
    yield assert_equal, str(ds), "snapshot_200"
    dso = [None, ("sphere", ("c", (0.1, 'unitary')))]
    dd = ds.all_data()
    yield assert_equal, dd["particle_position"].shape[0], 11907080
    yield assert_equal, dd["particle_position"].shape[1], 3
    tot = sum(dd[ptype, "particle_position"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    yield assert_equal, tot, 11907080
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        gdg, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(gdg, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2
