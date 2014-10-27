"""
Tipsy tests using the AGORA dataset




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.tipsy.api import TipsyDataset

_fields = (("deposit", "all_density"),
           ("deposit", "all_count"),
           ("deposit", "DarkMatter_density"),
)

pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
@requires_ds(pkdgrav, big_data = True, file_check = True)
def test_pkdgrav():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(field_dtypes = {"Coordinates": "d"},
                  cosmology_parameters = cosmology_parameters,
                  unit_base = {'length': (1.0/60.0, "Mpccm/h")},
                  n_ref = 64)
    ds = data_dir_load(pkdgrav, TipsyDataset, (), kwargs)
    yield assert_equal, str(ds), "halo1e11_run1.00400"
    dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
    dd = ds.all_data()
    yield assert_equal, dd["Coordinates"].shape, (26847360, 3)
    tot = sum(dd[ptype,"Coordinates"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    yield assert_equal, tot, 26847360
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        ds, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(ds, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2

gasoline = "agora_1e11.00400/agora_1e11.00400"
@requires_ds(gasoline, big_data = True, file_check = True)
def test_gasoline():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(cosmology_parameters = cosmology_parameters,
                  unit_base = {'length': (1.0/60.0, "Mpccm/h")},
                  n_ref = 64)
    ds = data_dir_load(gasoline, TipsyDataset, (), kwargs)
    yield assert_equal, str(ds), "agora_1e11.00400"
    dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
    dd = ds.all_data()
    yield assert_equal, dd["Coordinates"].shape, (10550576, 3)
    tot = sum(dd[ptype,"Coordinates"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    yield assert_equal, tot, 10550576
    for dobj_name in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "density"]:
                    yield PixelizedProjectionValuesTest(
                        ds, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(ds, field, dobj_name)
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2


@requires_file(pkdgrav)
def test_TipsyDataset():
    assert isinstance(data_dir_load(pkdgrav), TipsyDataset)
