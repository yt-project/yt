import numpy as np # NOQA
import pytest

import yt # NOQA
from yt.testing import requires_file
from yt.frontends.amrvac.api import AMRVACDataset, AMRVACGrid
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


blastwave_spherical_2D = "amrvac/bw_2d0000.dat"
khi_cartesian_2D = "amrvac/kh_2D0000.dat"
khi_cartesian_3D = "amrvac/kh_3D0000.dat"
jet_cylindrical_25D = "amrvac/Jet0003.dat"
riemann_cartesian_175D = "amrvac/R_1d0005.dat"
blastwave_cartesian_3D = "amrvac/bw_3d0000.dat"
blastwave_polar_2D = "amrvac/bw_polar_2D0000.dat"
blastwave_cylindrical_3D = "amrvac/bw_cylindrical_3D0000.dat"


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAMRVac(fw.AnswerTest):
    @requires_file(khi_cartesian_2D)
    def test_AMRVACDataset(self):
        assert isinstance(utils.data_dir_load(khi_cartesian_2D), AMRVACDataset)

    @utils.requires_ds(blastwave_cartesian_3D)
    def test_domain_size(self):
        #"Check for correct box size, see bw_3d.par"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        for lb in ds.domain_left_edge:
            assert int(lb) == 0
        for rb in ds.domain_right_edge:
            assert int(rb) == 2
        for w in ds.domain_width:
            assert int(w) == 2

    @requires_file(blastwave_cartesian_3D)
    def test_grid_attributes(self):
        #"Check various grid attributes"
        ds = utils.data_dir_load(blastwave_cartesian_3D)
        grids = ds.index.grids
        assert ds.index.max_level == 2
        for g in grids:
            assert isinstance(g, AMRVACGrid)
            assert isinstance(g.LeftEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.RightEdge, yt.units.yt_array.YTArray)
            assert isinstance(g.ActiveDimensions, np.ndarray)
            assert isinstance(g.Level, (np.int32, np.int64, int))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_polar_2D)
    def test_bw_polar_2d(self, a, d, w, f, ds_bw_polar_2D):
        self.hashes.update(self.small_patch_amr(ds_bw_polar_2d, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_cartesian_3D)
    def test_blastwave_cartesian_3D(self, a, d, w, f, ds_blastwave_cartesian_3D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_cartesian_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_spherical_2D)
    def test_blastwave_spherical_2D(self, a, d, w, f, ds_blastwave_spherical_2D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_spherical_2D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(blastwave_cylindrical_3D)
    def test_blastwave_cylindrical_3D(self, a, d, w, f, ds_blastwave_cylindrical_3D):
        self.hashes.update(self.small_patch_amr(ds_blastwave_cylindrical_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(khi_cartesian_2D)
    def test_khi_cartesian_2D(self, a, d, w, f, ds_khi_cartesian_2D):
        self.hashes.update(self.small_patch_amr(ds_khi_cartesian_2D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(khi_cartesian_3D)
    def test_khi_cartesian_3D(self, a, d, w, f, ds_khi_cartesian_3D):
        self.hashes.update(self.small_patch_amr(ds_khi_cartesian_3D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(jet_cylindrical_25D)
    def test_jet_cylindrical_25D(self, a, d, w, f, ds_jet_cylindrical_25D):
        self.hashes.update(self.small_patch_amr(ds_jet_cylindrical_25D, f, w, a, d))

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(riemann_cartesian_175D)
    def test_riemann_cartesian_175D(self, a, d, w, f, ds_riemann_cartesian_175D):
        self.hashes.update(self.small_patch_amr(ds_riemann_cartesian_175D, f, w, a, d))
