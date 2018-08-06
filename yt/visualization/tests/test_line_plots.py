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
from yt.testing import fake_random_ds
import os
import tempfile
import shutil

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, plot, test_prefix, decimals=12):
    def image_from_plot(filename_prefix):
        return plot.save(filename_prefix)

    image_from_plot.__name__ = "line_{}".format(test_prefix)
    test = GenericImageTest(ds, image_from_plot, decimals)
    test.prefix = test_prefix
    return test

tri2 = "SecondOrderTris/RZ_p_no_parts_do_nothing_bcs_cone_out.e"

@requires_ds(tri2)
def test_line_plot():
    ds = data_dir_load(tri2, kwargs={'step':-1})
    fields = [field for field in ds.field_list if field[0] == 'all']
    plot = yt.LinePlot(ds, fields, (0, 0, 0), (1, 1, 0), 1000)
    yield compare(ds, plot, "answers_line_plot")

@requires_ds(tri2)
def test_multi_line_plot():
    ds = data_dir_load(tri2, kwargs={'step':-1})
    fields = [field for field in ds.field_list if field[0] == 'all']
    lines = []
    lines.append(yt.LineBuffer(ds, [0.25, 0, 0], [0.25, 1, 0], 100, label='x = 0.25'))
    lines.append(yt.LineBuffer(ds, [0.5, 0, 0], [0.5, 1, 0], 100, label='x = 0.5'))
    plot = yt.LinePlot.from_lines(ds, fields, lines)
    yield compare(ds, plot, "answers_multi_line_plot")

def test_line_plot_methods():
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = fake_random_ds(32)

    plot = yt.LinePlot(ds, 'density', [0, 0, 0], [1, 1, 1], 512)
    plot.annotate_legend('density')
    plot.set_x_unit('cm')
    plot.set_unit('density', 'kg/cm**3')
    plot.save()

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)

def test_line_buffer():
    ds = fake_random_ds(32)
    lb = yt.LineBuffer(ds, (0, 0, 0), (1, 1, 1), 512, label='diag')
    lb['density']
    lb['density'] = 0
    lb['velocity_x']
    lb.keys()
    del lb['velocity_x']
