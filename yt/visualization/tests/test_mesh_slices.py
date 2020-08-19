"""
Tests for making unstructured mesh slices

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
import os
import tempfile

import numpy as np
import pytest

import yt
from yt.testing import small_fake_hexahedral_ds
from yt.utilities.answer_testing.answer_tests import generic_image_test
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.mesh_triangulation import triangulate_indices

# These tests randomly fail, so we skip them for now
pytest.skip("Mesh slice tests randomly fail. Skipping.", allow_module_level=True)


def slice_image(ds, field, idir):
    tmpfd, tmpfname = tempfile.mkstemp(suffix=".png")
    os.close(tmpfd)
    sl = yt.SlicePlot(ds, idir, field)
    sl.set_log("all", False)
    sl.save(tmpfname)
    return tmpfname


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir", "hashing")
class TestMesh:
    def test_mesh_slices_amr(self, ds_amr, field):
        img_fname = slice_image(ds_amr, field, 0)
        gi = generic_image_test(img_fname)
        self.hashes.update({"generic_image": gi})

    def test_mesh_slices_tetrahedral(self, ds_tetra, field, idir):
        mesh = ds_tetra.index.meshes[0]
        ad = ds_tetra.all_data()
        img_fname = slice_image(ds_tetra, field, idir)
        gi = generic_image_test(img_fname)
        self.hashes.update({"generic_image": gi})
        sl_obj = ds_tetra.slice(idir, ds_tetra.domain_center[idir])
        assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
        assert sl_obj[field].shape[0] < ad[field].shape[0]

    def test_mesh_slices_hexahedral(self, ds_hex, field, idir):
        # hexahedral ds
        ad = ds_hex.all_data()
        mesh = ds_hex.index.meshes[0]
        img_fname = slice_image(ds_hex, field, idir)
        gi = generic_image_test(img_fname)
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
