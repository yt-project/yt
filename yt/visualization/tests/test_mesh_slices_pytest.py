import numpy as np
import pytest

import yt
from yt.testing import small_fake_hexahedral_ds
from yt.utilities.answer_testing.answer_tests import generic_image
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.mesh_triangulation import triangulate_indices


def compare(ds, field, idir):
    def slice_image(im_name):
        sl = yt.SlicePlot(ds, idir, field)
        sl.set_log("all", False)
        image_file = sl.save(im_name)
        return image_file

    gi = generic_image(slice_image)
    # generic_image returns a list. In this case, there's only one entry,
    # which is a np array with the data we want
    return gi[0]


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir", "hashing")
class TestMesh:
    answer_file = None
    saved_hashes = None

    def test_mesh_slices_amr(self, ds_amr, field):
        gi = compare(ds_amr, field, 0)
        self.hashes.update({"generic_image": gi})

    def test_mesh_slices_tetrahedral(self, ds_tetra, field, idir):
        mesh = ds_tetra.index.meshes[0]
        ad = ds_tetra.all_data()
        gi = compare(ds_tetra, field, idir)
        self.hashes.update({"generic_image": gi})
        sl_obj = ds_tetra.slice(idir, ds_tetra.domain_center[idir])
        assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
        assert sl_obj[field].shape[0] < ad[field].shape[0]

    def test_mesh_slices_hexahedral(self, ds_hex, field, idir):
        # hexahedral ds
        ad = ds_hex.all_data()
        mesh = ds_hex.index.meshes[0]
        gi = compare(ds_hex, field, idir)
        self.hashes.update({"generic_image": gi})
        sl_obj = ds_hex.slice(idir, ds_hex.domain_center[idir])
        assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
        assert sl_obj[field].shape[0] < ad[field].shape[0]


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
