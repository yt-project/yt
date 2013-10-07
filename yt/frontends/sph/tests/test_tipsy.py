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
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    PixelizedProjectionValuesTest, \
    FieldValuesTest, \
    create_obj
from yt.frontends.sph.api import TipsyStaticOutput

_fields = (("deposit", "all_density"),
           ("deposit", "all_count"),
           ("deposit", "DarkMatter_density"),
)

pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
@requires_pf(pkdgrav, file_check = True)
def test_pkdgrav():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(endian="<",
                  field_dtypes = {"Coordinates": "d"},
                  cosmology_parameters = cosmology_parameters,
                  unit_base = {'mpchcm': 1.0/60.0},
                  n_ref = 64)
    pf = data_dir_load(pkdgrav, TipsyStaticOutput, (), kwargs)
    yield assert_equal, str(pf), "halo1e11_run1.00400"
    dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
    dd = pf.h.all_data()
    yield assert_equal, dd["Coordinates"].shape, (26847360, 3)
    tot = sum(dd[ptype,"Coordinates"].shape[0]
              for ptype in pf.particle_types if ptype != "all")
    yield assert_equal, tot, 26847360
    for ds in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        pf, axis, field, weight_field,
                        ds)
            yield FieldValuesTest(pf, field, ds)
        dobj = create_obj(pf, ds)
        s1 = dobj["Ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2

gasoline = "agora_1e11.00400/agora_1e11.00400"
@requires_pf(gasoline, file_check = True)
def test_gasoline():
    cosmology_parameters = dict(current_redshift = 0.0,
                                omega_lambda = 0.728,
                                omega_matter = 0.272,
                                hubble_constant = 0.702)
    kwargs = dict(cosmology_parameters = cosmology_parameters,
                  unit_base = {'mpchcm': 1.0/60.0},
                  n_ref = 64)
    pf = data_dir_load(gasoline, TipsyStaticOutput, (), kwargs)
    yield assert_equal, str(pf), "agora_1e11.00400"
    dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
    dd = pf.h.all_data()
    yield assert_equal, dd["Coordinates"].shape, (26847360, 3)
    tot = sum(dd[ptype,"Coordinates"].shape[0]
              for ptype in pf.particle_types if ptype != "all")
    yield assert_equal, tot, 26847360
    for ds in dso:
        for field in _fields:
            for axis in [0, 1, 2]:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        pf, axis, field, weight_field,
                        ds)
            yield FieldValuesTest(pf, field, ds)
        dobj = create_obj(pf, ds)
        s1 = dobj["Ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        yield assert_equal, s1, s2
