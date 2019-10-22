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

from yt.testing import \
    assert_allclose_units, \
    assert_equal, \
    disable_dataset_cache, \
    requires_file
from yt.frontends.athena.api import AthenaDataset
from yt.convenience import load
import yt.units as u
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
cloud = "ShockCloud/id0/Cloud.0050.vtk"
blast = "MHDBlast/id0/Blast.0100.vtk"
stripping = "RamPressureStripping/id0/rps.0062.vtk"
sloshing = "MHDSloshing/virgo_low_res.0054.vtk"


# Test data params
uo_stripping = {"time_unit":3.086e14,
                "length_unit":8.0236e22,
                "mass_unit":9.999e-30*8.0236e22**3}
uo_blast = {
    'length_unit': (1.0, 'pc'),
    'mass_unit': (2.38858753789e-24, 'g/cm**3*pc**3'),
    'time_unit': (1.0, 's*pc/km'),
}
uo_sloshing = {"length_unit": (1.0,"Mpc"),
               "time_unit": (1.0,"Myr"),
               "mass_unit": (1.0e14,"Msun")}


#============================================
#                 TestAthena
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('answer_file')
class TestAthena(fw.AnswerTest):
    #-----
    # test_cloud
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(cloud)
    def test_cloud(self, f, a, d, w, ds_cloud):
        ds = ds_cloud
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds, f, w, a, d))

    #-----
    # test_blast
    #-----
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blast)
    def test_blast(self, f, a, d, w, ds_blast):
        ds = ds_blast
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds, f, w, a, d))

    #-----
    # test_blast_override
    #-----
    @requires_file(blast)
    def test_blast_override(self):
        # verify that overriding units causes derived unit values to be updated.
        # see issue #1259
        ds = load(blast, units_override=uo_blast)
        assert_equal(float(ds.magnetic_unit.in_units('gauss')), 5.478674679698131e-07)

    #-----
    # test_stripping
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(stripping)
    def test_stripping(self, f, a, d, w, ds_stripping):
        # Run the small_patch_amr test suite
        self.hashes.update(self.small_patch_amr(ds, f, w, a, d))

    #-----
    # test_nprocs
    #-----
    @requires_file(sloshing)
    @disable_dataset_cache
    def test_nprocs(self):
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
        assert_equal(sp1.quantities.extrema("pressure"),
                     sp2.quantities.extrema("pressure"))
        assert_allclose_units(sp1.quantities.total_quantity("pressure"),
                              sp2.quantities.total_quantity("pressure"))
        for ax in "xyz":
            assert_equal(sp1.quantities.extrema("velocity_%s" % ax),
                         sp2.quantities.extrema("velocity_%s" % ax))
        assert_allclose_units(sp1.quantities.bulk_velocity(),
                              sp2.quantities.bulk_velocity())
        assert_equal(prj1["density"], prj2["density"])

    #-----
    # test_AthenaDataset
    #-----
    @requires_file(cloud)
    def test_AthenaDataset(self, ds_cloud):
        assert isinstance(ds_cloud, AthenaDataset)
