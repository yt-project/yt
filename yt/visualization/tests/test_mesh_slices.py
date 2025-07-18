import numpy as np
import pytest

import yt
from yt.testing import fake_hexahedral_ds, fake_tetrahedral_ds, small_fake_hexahedral_ds
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.mesh_triangulation import triangulate_indices


def setup_module():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class TestTetrahedral:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_tetrahedral_ds()
        cls.mesh = cls.ds.index.meshes[0]
        cls.ad = cls.ds.all_data()

        cls.slices = [
            cls.ds.slice(idir, cls.ds.domain_center[idir]) for idir in range(3)
        ]
        cls.sps = [yt.SlicePlot(cls.ds, idir, cls.ds.field_list) for idir in range(3)]
        for sp in cls.sps:
            sp.annotate_mesh_lines()
            sp.set_log("all", False)
            sp.render()

    @pytest.mark.parametrize("idir", range(3))
    def test_mesh_selection(self, idir):
        sl_obj = self.slices[idir]
        for field in self.ds.field_list:
            assert sl_obj[field].shape[0] == self.mesh.count(sl_obj.selector)
            assert sl_obj[field].shape[0] < self.ad[field].shape[0]

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_tetraheadral_all_elem(self, ax):
        return self.sps[ax].plots["all", "elem"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_tetraheadral_all_test(self, ax):
        return self.sps[ax].plots["all", "test"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_tetraheadral_connect1_elem(self, ax):
        return self.sps[ax].plots["connect1", "elem"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_tetraheadral_connect1_test(self, ax):
        return self.sps[ax].plots["connect1", "test"].figure


class TestHexahedral:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_hexahedral_ds()
        cls.mesh = cls.ds.index.meshes[0]
        cls.ad = cls.ds.all_data()

        cls.slices = [
            cls.ds.slice(idir, cls.ds.domain_center[idir]) for idir in range(3)
        ]
        cls.sps = [yt.SlicePlot(cls.ds, idir, cls.ds.field_list) for idir in range(3)]
        for sp in cls.sps:
            sp.annotate_mesh_lines()
            sp.set_log("all", False)
            sp.render()

    @pytest.mark.parametrize("idir", range(3))
    def test_mesh_selection(self, idir):
        sl_obj = self.slices[idir]
        for field in self.ds.field_list:
            assert sl_obj[field].shape[0] == self.mesh.count(sl_obj.selector)
            assert sl_obj[field].shape[0] < self.ad[field].shape[0]

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_hexaheadral_all_elem(self, ax):
        return self.sps[ax].plots["all", "elem"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_hexaheadral_all_test(self, ax):
        return self.sps[ax].plots["all", "test"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_hexaheadral_connect1_elem(self, ax):
        return self.sps[ax].plots["connect1", "elem"].figure

    @pytest.mark.parametrize("ax", range(3))
    @pytest.mark.mpl_image_compare
    def test_mesh_slice_hexaheadral_connect1_test(self, ax):
        return self.sps[ax].plots["connect1", "test"].figure


def test_perfect_element_intersection():
    # This test tests mesh line annotation where a z=0 slice
    # perfectly intersects the top of a hexahedral element with node
    # z-coordinates containing both -0 and +0. Before
    # https://github.com/yt-project/yt/pull/1437 this test falsely
    # yielded three annotation lines, whereas the correct result is four
    # corresponding to the four edges of the top hex face.

    ds = small_fake_hexahedral_ds()
    indices = ds.index.meshes[0].connectivity_indices
    coords = ds.index.meshes[0].connectivity_coords
    tri_indices = triangulate_indices(indices)
    tri_coords = coords[tri_indices]
    lines = triangle_plane_intersect(2, 0, tri_coords)
    non_zero_lines = 0
    for i in range(lines.shape[0]):
        norm = np.linalg.norm(lines[i][0] - lines[i][1])
        if norm > 1e-8:
            non_zero_lines += 1
    np.testing.assert_equal(non_zero_lines, 4)
