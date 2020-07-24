"""
Title: test_athena.py
Purpose: Athena frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.athena.api import AthenaDataset
from yt.convenience import load
from yt.testing import (
    assert_allclose_units,
    assert_equal,
    disable_dataset_cache,
    requires_file,
)
import yt.units as u
from yt.utilities.answer_testing.answer_tests import small_patch_amr
from yt.utilities.answer_testing.utils import requires_ds


# Test data
cloud = "ShockCloud/id0/Cloud.0050.vtk"
blast = "MHDBlast/id0/Blast.0100.vtk"
stripping = "RamPressureStripping/id0/rps.0062.vtk"
sloshing = "MHDSloshing/virgo_low_res.0054.vtk"


# Test data params
uo_stripping = {
    "time_unit": 3.086e14,
    "length_unit": 8.0236e22,
    "mass_unit": 9.999e-30 * 8.0236e22 ** 3,
}
uo_blast = {
    "length_unit": (1.0, "pc"),
    "mass_unit": (2.38858753789e-24, "g/cm**3*pc**3"),
    "time_unit": (1.0, "s*pc/km"),
}
uo_sloshing = {
    "length_unit": (1.0, "Mpc"),
    "time_unit": (1.0, "Myr"),
    "mass_unit": (1.0e14, "Msun"),
}


@pytest.mark.answer_test
@pytest.mark.usefixtures("answer_file")
class TestAthena:
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [cloud], indirect=True)
    def test_cloud(self, f, a, d, w, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [blast], indirect=True)
    def test_blast(self, f, a, d, w, ds):
        self.hashes.update(small_patch_amr(ds, f, w, a, d))

    @requires_file(blast)
    def test_blast_override(self):
        # verify that overriding units causes derived unit values to be updated.
        # see issue #1259
        ds = load(blast, units_override=uo_blast)
        assert_equal(float(ds.magnetic_unit.in_units("gauss")), 5.478674679698131e-07)

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @requires_ds(stripping)
    def test_stripping(self, f, a, d, w, ds_stripping):
        self.hashes.update(small_patch_amr(ds_stripping, f, w, a, d))

    @requires_file(sloshing)
    @disable_dataset_cache
    def test_nprocs(self):
        ds1 = load(sloshing, units_override=uo_sloshing)
        sp1 = ds1.sphere("c", (100.0, "kpc"))
        prj1 = ds1.proj("density", 0)
        ds2 = load(sloshing, units_override=uo_sloshing, nprocs=8)
        sp2 = ds2.sphere("c", (100.0, "kpc"))
        prj2 = ds1.proj("density", 0)
        ds3 = load(sloshing, parameters=uo_sloshing)
        assert_equal(ds3.length_unit, u.Mpc)
        assert_equal(ds3.time_unit, u.Myr)
        assert_equal(ds3.mass_unit, 1e14 * u.Msun)
        assert_equal(
            sp1.quantities.extrema("pressure"), sp2.quantities.extrema("pressure")
        )
        assert_allclose_units(
            sp1.quantities.total_quantity("pressure"),
            sp2.quantities.total_quantity("pressure"),
        )
        for ax in "xyz":
            assert_equal(
                sp1.quantities.extrema("velocity_%s" % ax),
                sp2.quantities.extrema("velocity_%s" % ax),
            )
        assert_allclose_units(
            sp1.quantities.bulk_velocity(), sp2.quantities.bulk_velocity()
        )
        assert_equal(prj1["density"], prj2["density"])

    @pytest.mark.parametrize("ds", [cloud], indirect=True)
    def test_AthenaDataset(self, ds):
        assert isinstance(ds, AthenaDataset)
