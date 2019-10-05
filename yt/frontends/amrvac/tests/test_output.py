import numpy as np

import yt
from yt.utilities.parameter_file_storage import output_type_registry
from yt.utilities.answer_testing.framework import requires_ds, data_dir_load

blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"
khi_cartesian_2D = "amrvac/kh_2D0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"

@requires_ds(blastwave_cartesian_3D)
def test_domain_size(self):
    """Check for correct box size, see bw_3d.par"""
    ds = data_dir_load(blastwave_cartesian_3D)
    for lb in ds.domain_left_edge:
        assert int(lb) == 0
    for rb in ds.domain_right_edge:
        assert int(rb) == 2
    for w in ds.domain_width:
        assert int(w) == 2

@requires_ds(blastwave_cartesian_3D)
def test_grid_attributes(self):
    """Check various grid attributes"""
    ds = data_dir_load(blastwave_cartesian_3D)
    grids = ds.index.grids
    assert ds.index.max_level == 2
    for g in grids:
        assert isinstance(g, yt.frontends.amrvac.AMRVACGrid)
        assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
        assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
        assert isinstance(g.ActiveDimensions, np.ndarray)
        assert isinstance(g.Level, (np.int32, np.int64, int))

