"""
Tests for making line plots

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import yt
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    GenericImageTest

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, fields, point1, point2, resolution, test_prefix, decimals=12):
    def line_plot(filename_prefix):
        ln = yt.LinePlot(ds, fields, point1, point2, resolution)
        image_file = ln.save(filename_prefix)
        return image_file

    line_plot.__name__ = "line_{}".format(test_prefix)
    test = GenericImageTest(ds, line_plot, decimals)
    test.prefix = test_prefix
    return test

tri2 = "SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e"
iso_galaxy = "IsolatedGalaxy/galaxy0030/galaxy0030"

@requires_ds(tri2)
def test_line_plot():
    ds = data_dir_load(tri2, kwargs={'step':-1})
    fields = [field for field in ds.field_list if field[0] == 'all']
    yield compare(ds, fields, (0, 0, 0), (1, 1, 0), 1000, "answers_line_plot")

@requires_ds(iso_galaxy)
def test_line_plot_methods():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = data_dir_load(iso_galaxy)

    plot = yt.LinePlot(ds, 'density', [0, 0, 0], [1, 1, 1], 512)
    plot.add_legend('density')
    plot.set_x_unit('cm')
    plot.set_unit('density', 'kg/cm**3')
    plot.save()

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
