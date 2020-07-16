from collections import OrderedDict

import numpy as np

from yt.frontends.flash.api import FLASHDataset, FLASHParticleDataset
from yt.testing import (
    ParticleSelectionComparison,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    nbody_answer,
    requires_ds,
    small_patch_amr,
)

_fields = ("temperature", "density", "velocity_magnitude")

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"


@requires_ds(sloshing, big_data=True)
def test_sloshing():
    ds = data_dir_load(sloshing)
    assert_equal(str(ds), "sloshing_low_res_hdf5_plt_cnt_0300")
    for test in small_patch_amr(ds, _fields):
        test_sloshing.__name__ = test.description
        yield test


_fields_2d = ("temperature", "density")

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
    dd["density"]


@requires_ds(fid_1to3_b1, big_data=True)
def test_fid_1to3_b1():
    ds = data_dir_load(fid_1to3_b1)
    for test in nbody_answer(
        ds, "fiducial_1to3_b1_hdf5_part_0080", 6684119, fid_1to3_b1_fields
    ):
        test_fid_1to3_b1.__name__ = test.description
        yield test
