import numpy as np
from nose.plugins.attrib import attr

import yt
from yt.testing import (
    ANSWER_TEST_TAG,
    fake_amr_ds,
    fake_hexahedral_ds,
    fake_tetrahedral_ds,
    small_fake_hexahedral_ds,
)
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.mesh_triangulation import triangulate_indices


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def compare(ds, field, idir, test_prefix, test_name, decimals=12, annotate=False):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, idir, field)
        if annotate:
            sl.annotate_mesh_lines()
        sl.set_log("all", False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = f"slice_{test_prefix}"
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test


@attr(ANSWER_TEST_TAG)
def test_mesh_slices_amr():
    ds = fake_amr_ds()
    for field in ds.field_list:
        prefix = f"{field[0]}_{field[1]}_{0}"
        yield compare(ds, field, 0, test_prefix=prefix, test_name="mesh_slices_amr")


@attr(ANSWER_TEST_TAG)
def test_mesh_slices_tetrahedral():
    ds = fake_tetrahedral_ds()

    mesh = ds.index.meshes[0]
    ad = ds.all_data()

    for field in ds.field_list:
        for idir in [0, 1, 2]:
            prefix = f"{field[0]}_{field[1]}_{idir}"
            yield compare(
                ds,
                field,
                idir,
                test_prefix=prefix,
                test_name="mesh_slices_tetrahedral",
                annotate=True,
            )

            sl_obj = ds.slice(idir, ds.domain_center[idir])
            assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
            assert sl_obj[field].shape[0] < ad[field].shape[0]


@attr(ANSWER_TEST_TAG)
def test_mesh_slices_hexahedral():
    # hexahedral ds
    ds = fake_hexahedral_ds()
    ad = ds.all_data()
    mesh = ds.index.meshes[0]

    for field in ds.field_list:
        for idir in [0, 1, 2]:
            prefix = f"{field[0]}_{field[1]}_{idir}"
            yield compare(
                ds,
                field,
                idir,
                test_prefix=prefix,
                test_name="mesh_slices_hexahedral",
                annotate=True,
            )

            sl_obj = ds.slice(idir, ds.domain_center[idir])
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
