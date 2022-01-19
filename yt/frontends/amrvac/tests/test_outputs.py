import numpy as np

import yt  # NOQA
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
from yt.testing import requires_file
from yt.units import YTArray
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2d0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"
rmi_cartesian_dust_2D = "amrvac/Richtmyer_Meshkov_dust_2D/RM2D_dust_Kwok0000.dat"


def _get_fields_to_check(ds):
    fields = ["density", "velocity_magnitude"]
    raw_fields_labels = [fname for ftype, fname in ds.field_list]
    if "b1" in raw_fields_labels:
        fields.append("magnetic_energy_density")
    if "e" in raw_fields_labels:
        fields.append("energy_density")
    if "rhod1" in raw_fields_labels:
        fields.append("total_dust_density")
        # note : not hitting dust velocity fields
    return fields


@requires_file(khi_cartesian_2D)
def test_AMRVACDataset():
    assert isinstance(data_dir_load(khi_cartesian_2D), AMRVACDataset)


@requires_ds(blastwave_cartesian_3D)
def test_domain_size():
    # "Check for correct box size, see bw_3d.par"
    ds = data_dir_load(blastwave_cartesian_3D)
    for lb in ds.domain_left_edge:
        assert int(lb) == 0
    for rb in ds.domain_right_edge:
        assert int(rb) == 2
    for w in ds.domain_width:
        assert int(w) == 2


@requires_file(blastwave_cartesian_3D)
def test_grid_attributes():
    # "Check various grid attributes"
    ds = data_dir_load(blastwave_cartesian_3D)
    grids = ds.index.grids
    assert ds.index.max_level == 2
    for g in grids:
        assert isinstance(g, AMRVACGrid)
        assert isinstance(g.LeftEdge, YTArray)
        assert isinstance(g.RightEdge, YTArray)
        assert isinstance(g.ActiveDimensions, np.ndarray)
        assert isinstance(g.Level, (np.int32, np.int64, int))


@requires_ds(blastwave_polar_2D)
def test_bw_polar_2d():
    ds = data_dir_load(blastwave_polar_2D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_bw_polar_2d.__name__ = test.description
        yield test


@requires_ds(blastwave_cartesian_3D)
def test_blastwave_cartesian_3D():
    ds = data_dir_load(blastwave_cartesian_3D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_blastwave_cartesian_3D.__name__ = test.description
        yield test


@requires_ds(blastwave_spherical_2D)
def test_blastwave_spherical_2D():
    ds = data_dir_load(blastwave_spherical_2D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_blastwave_spherical_2D.__name__ = test.description
        yield test


@requires_ds(blastwave_cylindrical_3D)
def test_blastwave_cylindrical_3D():
    ds = data_dir_load(blastwave_cylindrical_3D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_blastwave_cylindrical_3D.__name__ = test.description
        yield test


@requires_ds(khi_cartesian_2D)
def test_khi_cartesian_2D():
    ds = data_dir_load(khi_cartesian_2D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_khi_cartesian_2D.__name__ = test.description
        yield test


@requires_ds(khi_cartesian_3D)
def test_khi_cartesian_3D():
    ds = data_dir_load(khi_cartesian_3D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_khi_cartesian_3D.__name__ = test.description
        yield test


@requires_ds(jet_cylindrical_25D)
def test_jet_cylindrical_25D():
    ds = data_dir_load(jet_cylindrical_25D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_jet_cylindrical_25D.__name__ = test.description
        yield test


@requires_ds(riemann_cartesian_175D)
def test_riemann_cartesian_175D():
    ds = data_dir_load(riemann_cartesian_175D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_riemann_cartesian_175D.__name__ = test.description
        yield test


@requires_ds(rmi_cartesian_dust_2D)
def test_rmi_cartesian_dust_2D():
    # dataset with dust fields
    ds = data_dir_load(rmi_cartesian_dust_2D)
    for test in small_patch_amr(ds, _get_fields_to_check(ds)):
        test_rmi_cartesian_dust_2D.__name__ = test.description
        yield test
