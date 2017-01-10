"""
Athena frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    assert_equal, \
    requires_file, \
    assert_allclose_units
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load
from yt.frontends.athena.api import AthenaDataset
from yt.config import ytcfg
from yt.convenience import load

import yt.units as u

_fields_cloud = ("scalar[0]", "density", "total_energy")

cloud = "ShockCloud/id0/Cloud.0050.vtk"
@requires_ds(cloud)
def test_cloud():
    ds = data_dir_load(cloud)
    yield assert_equal, str(ds), "Cloud.0050"
    for test in small_patch_amr(ds, _fields_cloud):
        test_cloud.__name__ = test.description
        yield test

_fields_blast = ("temperature", "density", "velocity_magnitude")

blast = "MHDBlast/id0/Blast.0100.vtk"
@requires_ds(blast)
def test_blast():
    ds = data_dir_load(blast)
    yield assert_equal, str(ds), "Blast.0100"
    for test in small_patch_amr(ds, _fields_blast):
        test_blast.__name__ = test.description
        yield test

uo_blast = {
    'length_unit': (1.0, 'pc'),
    'mass_unit': (2.38858753789e-24, 'g/cm**3*pc**3'),
    'time_unit': (1.0, 's*pc/km'),
}

@requires_file(blast)
def test_blast_override():
    # verify that overriding units causes derived unit values to be updated.
    # see issue #1259
    ds = load(blast, units_override=uo_blast)
    assert_equal(float(ds.magnetic_unit.in_units('gauss')),
                 5.478674679698131e-07)

uo_stripping = {"time_unit":3.086e14,
                "length_unit":8.0236e22,
                "mass_unit":9.999e-30*8.0236e22**3}

_fields_stripping = ("temperature", "density", "specific_scalar[0]")

stripping = "RamPressureStripping/id0/rps.0062.vtk"
@requires_ds(stripping, big_data=True)
def test_stripping():
    ds = data_dir_load(stripping, kwargs={"units_override":uo_stripping})
    yield assert_equal, str(ds), "rps.0062"
    for test in small_patch_amr(ds, _fields_stripping):
        test_stripping.__name__ = test.description
        yield test

sloshing = "MHDSloshing/virgo_low_res.0054.vtk"

uo_sloshing = {"length_unit": (1.0,"Mpc"),
               "time_unit": (1.0,"Myr"),
               "mass_unit": (1.0e14,"Msun")}

@requires_file(sloshing)
def test_nprocs():
    ytcfg["yt","skip_dataset_cache"] = "True"

    ds1 = load(sloshing, units_override=uo_sloshing)
    sp1 = ds1.sphere("c", (100.,"kpc"))
    prj1 = ds1.proj("density",0)
    ds2 = load(sloshing, units_override=uo_sloshing, nprocs=8)
    sp2 = ds2.sphere("c", (100.,"kpc"))
    prj2 = ds1.proj("density",0)

    ds3 = load(sloshing, parameters=uo_sloshing)
    assert_equal(ds3.length_unit, u.Mpc)
    assert_equal(ds3.time_unit, u.Myr)
    assert_equal(ds3.mass_unit, 1e14*u.Msun)

    yield assert_equal, sp1.quantities.extrema("pressure"), sp2.quantities.extrema("pressure")
    yield assert_allclose_units, sp1.quantities.total_quantity("pressure"), sp2.quantities.total_quantity("pressure")
    for ax in "xyz":
        yield assert_equal, sp1.quantities.extrema("velocity_%s" % ax), sp2.quantities.extrema("velocity_%s" % ax)
    yield assert_allclose_units, sp1.quantities.bulk_velocity(), sp2.quantities.bulk_velocity()
    yield assert_equal, prj1["density"], prj2["density"]

    ytcfg["yt","skip_dataset_cache"] = "False"

@requires_file(cloud)
def test_AthenaDataset():
    assert isinstance(data_dir_load(cloud), AthenaDataset)
