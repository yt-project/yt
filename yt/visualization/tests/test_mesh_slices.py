"""
Tests for making unstructured mesh slices

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from collections import OrderedDict
import os
import tempfile

import numpy as np
import pytest

import yt
from yt.testing import fake_amr_ds, fake_tetrahedral_ds, \
    fake_hexahedral_ds, small_fake_hexahedral_ds
import yt.utilities.answer_testing.framework as fw 
from yt.utilities.answer_testing import utils
from yt.utilities.lib.geometry_utils import triangle_plane_intersect
from yt.utilities.lib.mesh_triangulation import triangulate_indices


# Answer file
answer_file = 'mesh_slices_answers.yaml'


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def slice_image(ds, field, idir):
    tmpfd, tmpfname = tempfile.mkstemp(suffix='.png')
    os.close(tmpfd)
    sl = yt.SlicePlot(ds, idir, field)
    sl.set_log('all', False)
    sl.save(tmpfname)
    return tmpfname 


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir')
class TestMesh(fw.AnswerTest):
    def test_mesh_slices_amr(self):
        ds = fake_amr_ds()
        hd = OrderedDict()
        hd['generic_image'] = OrderedDict()
        for field in ds.field_list:
            img_fname = slice_image(ds, field, 0)
            gi_hd = utils.generate_hash(self.generic_image_test(img_fname))
            hd['generic_image'][field] = gi_hd
        hd = {'mesh_slices_amr' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hd, self.answer_store)

    def test_mesh_slices_tetrahedral(self):
        ds = fake_tetrahedral_ds()
        mesh = ds.index.meshes[0]
        ad = ds.all_data()
        hd = OrderedDict()
        hd['generic_image'] = OrderedDict()
        for field in ds.field_list:
            hd['generic_image'][field] = OrderedDict()
            for idir in [0, 1, 2]:
                img_fname = slice_image(ds, field, idir)
                gi_hd = utils.generate_hash(self.generic_image_test(img_fname))
                hd['generic_image'][field][str(idir)] = gi_hd
                sl_obj = ds.slice(idir, ds.domain_center[idir])
                assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
                assert sl_obj[field].shape[0] < ad[field].shape[0]
        hd = {'mesh_slices_tetrahedral' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hd, self.answer_store)

    def test_mesh_slices_hexahedral(self):
        # hexahedral ds
        ds = fake_hexahedral_ds()
        ad = ds.all_data()
        mesh = ds.index.meshes[0]
        hd = OrderedDict()
        hd['generic_image'] = OrderedDict()
        for field in ds.field_list:
            hd['generic_image'][field] = OrderedDict()
            for idir in [0, 1, 2]:
                img_fname = slice_image(ds, field, idir)
                gi_hd = utils.generate_hash(self.generic_image_test(img_fname))
                hd['generic_image'][field][str(idir)] = gi_hd
                sl_obj = ds.slice(idir, ds.domain_center[idir])
                assert sl_obj[field].shape[0] == mesh.count(sl_obj.selector)
                assert sl_obj[field].shape[0] < ad[field].shape[0]
        hd = {'mesh_slices_hexahedral' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hd, self.answer_store)

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
