import numpy as np

import yt
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, small_patch_amr
from yt.testing import requires_file
from yt.frontends.amrvac.api import AMRVACDataset

blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2D0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"

# the following tests are not yet uploaded to yt-website ! (wip)
#blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
#blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
#blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"

@requires_file(khi_cartesian_2D)
def test_AMRVACDataset():
    assert isinstance(data_dir_load(khi_cartesian_2D), AMRVACDataset)


#@requires_ds(blastwave_cartesian_3D)
#def test_domain_size():
#    #"Check for correct box size, see bw_3d.par"
#    ds = data_dir_load(blastwave_cartesian_3D)
#    for lb in ds.domain_left_edge:
#        assert int(lb) == 0
#    for rb in ds.domain_right_edge:
#        assert int(rb) == 2
#    for w in ds.domain_width:
#        assert int(w) == 2

#@requires_file(blastwave_cartesian_3D)
#def test_grid_attributes(self):
#    #"Check various grid attributes"
#    ds = data_dir_load(blastwave_cartesian_3D)
#    grids = ds.index.grids
#    assert ds.index.max_level == 2
#    for g in grids:
#        assert isinstance(g, yt.frontends.amrvac.AMRVACGrid)
#        assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
#        assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
#        assert isinstance(g.ActiveDimensions, np.ndarray)
#        assert isinstance(g.Level, (np.int32, np.int64, int))

#@requires_ds(blastwave_polar_2D)
#def test_bw_polar_2d():
#    ds = data_dir_load(blastwave_polar_2D)
#    for test in small_patch_amr(ds, ds.field_list):
#        test_bw_polar_2d.__name__ = test.description
#        yield test

#@requires_ds(blastwave_cartesian_3D)
#def test_blastwave_cartesian_3D():
#    ds = data_dir_load(blastwave_cartesian_3D)
#    for test in small_patch_amr(ds, ds.field_list):
#        test_blastwave_cartesian_3D.__name__ = test.description
#        yield test

@requires_ds(blastwave_spherical_2D)
def test_blastwave_spherical_2D():
    ds = data_dir_load(blastwave_spherical_2D)
    for test in small_patch_amr(ds, ds.field_list):
        test_blastwave_spherical_2D.__name__ = test.description
        yield test

#@requires_ds(blastwave_cylindrical_3D)
#def test_blastwave_cylindrical_3D():
#    ds = data_dir_load(blastwave_cylindrical_3D)
#    for test in small_patch_amr(ds, ds.field_list):
#        test_blastwave_cylindrical_3D.__name__ = test.description
#        yield test

@requires_ds(khi_cartesian_2D)
def test_khi_cartesian_2D():
    ds = data_dir_load(khi_cartesian_2D)
    for test in small_patch_amr(ds, ds.field_list):
        test_khi_cartesian_2D.__name__ = test.description
        yield test

@requires_ds(khi_cartesian_3D)
def test_khi_cartesian_3D():
    ds = data_dir_load(khi_cartesian_3D)
    for test in small_patch_amr(ds, ds.field_list):
        test_khi_cartesian_3D.__name__ = test.description
        yield test

#@requires_ds(jet_cylindrical_25D)
#def test_jet_cylindrical_25D():
#    ds = data_dir_load(jet_cylindrical_25D)
#    for test in small_patch_amr(ds, ds.field_list):
#        test_jet_cylindrical_25D.__name__ = test.description
#        yield test

@requires_ds(riemann_cartesian_175D)
def test_riemann_cartesian_175D():
    ds = data_dir_load(riemann_cartesian_175D)
    for test in small_patch_amr(ds, ds.field_list):
        test_riemann_cartesian_175D.__name__ = test.description
        yield test
