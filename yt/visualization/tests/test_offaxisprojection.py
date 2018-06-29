"""
Test for off_axis_projection and write_projection



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os
import os.path
import tempfile
import shutil
from yt.testing import \
    fake_random_ds, assert_equal, expand_keywords, fake_octree_ds
from yt.mods import write_projection
from yt.visualization.volume_rendering.api import off_axis_projection
from yt.visualization.api import OffAxisProjectionPlot, \
    OffAxisSlicePlot

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def test_oap(tmpdir=True):
    """Tests functionality of off_axis_projection and write_projection."""
    # Perform I/O in safe place instead of yt main dir
    if tmpdir:
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)

    # args for off_axis_projection
    test_ds = fake_random_ds(64)
    c = test_ds.domain_center 
    norm = [0.5, 0.5, 0.5]
    W = test_ds.arr([0.5,0.5,1.0], 'unitary')
    N = 256
    field = ("gas","density")
    oap_args = [test_ds, c, norm, W, N, field]

    # kwargs for off_axis_projection
    oap_kwargs = {}
    oap_kwargs['weight'] = (None, 'cell_mass')
    oap_kwargs['no_ghost'] = (True, False)
    oap_kwargs['interpolated'] = (False,)
    oap_kwargs['north_vector'] = ((1,0,0), (0,0.5,1.0))
    oap_kwargs_list = expand_keywords(oap_kwargs)

    # args or write_projection
    fn = "test_%d.png"

    # kwargs for write_projection
    wp_kwargs = {}
    wp_kwargs['colorbar'] = (True, False)
    wp_kwargs['colorbar_label'] = ('test')
    wp_kwargs['title'] = ('test')
    wp_kwargs['limits'] = (None, (1e3, 1e5))
    wp_kwargs['take_log'] = (True, False)
    wp_kwargs['figsize'] = ((8,6), [1,1])
    wp_kwargs['dpi'] = (100, 50)
    wp_kwargs['cmap_name'] = ('arbre', 'kelp')
    wp_kwargs_list = expand_keywords(wp_kwargs)

    # test all off_axis_projection kwargs and write_projection kwargs
    # make sure they are able to be projected, then remove and try next
    # iteration
    for i, oap_kwargs in enumerate(oap_kwargs_list):
        image = off_axis_projection(*oap_args, **oap_kwargs)
        for wp_kwargs in wp_kwargs_list:
            write_projection(image, fn % i, **wp_kwargs)
            assert_equal(os.path.exists(fn % i), True)

    if tmpdir:
        os.chdir(curdir)
        # clean up
        shutil.rmtree(tmpdir)

def test_field_cut_off_axis_octree():
    ds = fake_octree_ds() 
    cut = ds.all_data().cut_region('obj["density"]>0.5') 
    p1 = OffAxisProjectionPlot(ds, [1, 0, 0], 'density')
    p2 = OffAxisProjectionPlot(ds, [1, 0, 0], 'density', data_source=cut)
    assert_equal(p2.frb["density"].min() == 0.0, True) # Lots of zeros
    assert_equal((p1.frb["density"] == p2.frb["density"]).all(), False)
    p3 = OffAxisSlicePlot(ds, [1, 0, 0], 'density')
    p4 = OffAxisSlicePlot(ds, [1, 0, 0], 'density', data_source=cut)
    assert_equal((p3.frb["density"] == p4.frb["density"]).all(), False)
    p4rho = p4.frb["density"]
    assert_equal(p4rho.min() == 0.0, True) # Lots of zeros
    assert_equal(p4rho[p4rho > 0.0].min() >= 0.5, True)

if __name__ == "__main__":
    for test in test_oap(tmpdir=False):
        pass
