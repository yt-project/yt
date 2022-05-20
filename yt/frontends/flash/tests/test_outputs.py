from collections import OrderedDict

import numpy as np

from yt.frontends.flash.api import FLASHDataset, FLASHParticleDataset
from yt.loaders import load
from yt.testing import (
    ParticleSelectionComparison,
    assert_allclose,
    assert_equal,
    disable_dataset_cache,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    nbody_answer,
    requires_ds,
    small_patch_amr,
)

_fields = (("gas", "temperature"), ("gas", "density"), ("gas", "velocity_magnitude"))

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"


@requires_ds(sloshing, big_data=True)
def test_sloshing():
    ds = data_dir_load(sloshing)
    assert_equal(str(ds), "sloshing_low_res_hdf5_plt_cnt_0300")
    for test in small_patch_amr(ds, _fields):
        test_sloshing.__name__ = test.description
        yield test


_fields_2d = (("gas", "temperature"), ("gas", "density"))

wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"


@requires_ds(wt)
def test_wind_tunnel():
    ds = data_dir_load(wt)
    assert_equal(str(ds), "windtunnel_4lev_hdf5_plt_cnt_0030")
    for test in small_patch_amr(ds, _fields_2d):
        test_wind_tunnel.__name__ = test.description
        yield test


@requires_file(wt)
def test_FLASHDataset():
    assert isinstance(data_dir_load(wt), FLASHDataset)


@requires_file(sloshing)
def test_units_override():
    units_override_check(sloshing)


@disable_dataset_cache
@requires_file(sloshing)
def test_magnetic_units():
    ds1 = load(sloshing)
    assert_allclose(ds1.magnetic_unit.value, np.sqrt(4.0 * np.pi))
    assert str(ds1.magnetic_unit.units) == "G"
    mag_unit1 = ds1.magnetic_unit.to("code_magnetic")
    assert_allclose(mag_unit1.value, 1.0)
    assert str(mag_unit1.units) == "code_magnetic"
    ds2 = load(sloshing, unit_system="mks")
    assert_allclose(ds2.magnetic_unit.value, np.sqrt(4.0 * np.pi) * 1.0e-4)
    assert str(ds2.magnetic_unit.units) == "T"
    mag_unit2 = ds2.magnetic_unit.to("code_magnetic")
    assert_allclose(mag_unit2.value, 1.0)
    assert str(mag_unit2.units) == "code_magnetic"


@requires_file(sloshing)
def test_mu():
    ds = data_dir_load(sloshing)
    sp = ds.sphere("c", (0.1, "unitary"))
    assert np.all(
        sp["gas", "mean_molecular_weight"] == ds.parameters["eos_singlespeciesa"]
    )


fid_1to3_b1 = "fiducial_1to3_b1/fiducial_1to3_b1_hdf5_part_0080"

fid_1to3_b1_fields = OrderedDict(
    [
        (("all", "particle_mass"), None),
        (("all", "particle_ones"), None),
        (("all", "particle_velocity_x"), ("all", "particle_mass")),
        (("all", "particle_velocity_y"), ("all", "particle_mass")),
        (("all", "particle_velocity_z"), ("all", "particle_mass")),
    ]
)


@requires_file(fid_1to3_b1)
def test_FLASHParticleDataset():
    assert isinstance(data_dir_load(fid_1to3_b1), FLASHParticleDataset)


@requires_file(fid_1to3_b1)
def test_FLASHParticleDataset_selection():
    ds = data_dir_load(fid_1to3_b1)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


dens_turb_mag = "DensTurbMag/DensTurbMag_hdf5_plt_cnt_0015"


@requires_file(dens_turb_mag)
def test_FLASH25_dataset():
    ds = data_dir_load(dens_turb_mag)
    assert_equal(ds.parameters["time"], 751000000000.0)
    assert_equal(ds.domain_dimensions, np.array([8, 8, 8]))
    assert_equal(ds.domain_left_edge, ds.arr([-2e18, -2e18, -2e18], "code_length"))

    assert_equal(ds.index.num_grids, 73)
    dd = ds.all_data()
    dd[("gas", "density")]


@requires_ds(fid_1to3_b1, big_data=True)
def test_fid_1to3_b1():
    ds = data_dir_load(fid_1to3_b1)
    for test in nbody_answer(
        ds, "fiducial_1to3_b1_hdf5_part_0080", 6684119, fid_1to3_b1_fields
    ):
        test_fid_1to3_b1.__name__ = test.description
        yield test


loc_bub_dust = "LocBub_dust/LocBub_dust_hdf5_plt_cnt_0220"


@requires_file(loc_bub_dust)
def test_blockless_particles():
    ds = data_dir_load(loc_bub_dust)
    dd = ds.all_data()
    pos = dd["all", "particle_position"]
    assert_equal(pos.shape, (2239, 3))
