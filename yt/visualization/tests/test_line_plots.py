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
from nose.plugins.attrib import attr
from nose.tools import assert_raises

import yt
from yt.testing import assert_equal, fake_random_ds, ANSWER_TEST_TAG
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.line_plot import _validate_point

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, plot, test_prefix, test_name, decimals=12):
    def image_from_plot(filename_prefix):
        return plot.save(filename_prefix)

    image_from_plot.__name__ = "line_{}".format(test_prefix)
    test = GenericImageTest(ds, image_from_plot, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test

@attr(ANSWER_TEST_TAG)
def test_line_plot():
    ds = fake_random_ds(4)
    fields = [field for field in ds.field_list if field[0] == 'stream']
    field_labels = {f:f[1] for f in fields}
    plot = yt.LinePlot(ds, fields, (0, 0, 0), (1, 1, 0), 1000,
                       field_labels=field_labels)
    plot.annotate_legend(fields[0])
    plot.annotate_legend(fields[1])
    plot.set_x_unit('cm')
    plot.set_unit(fields[0], 'kg/cm**3')
    plot.annotate_title(fields[0], "Density Plot")
    yield compare(ds, plot, test_prefix="answers_line_plot",
                  test_name="answers_line_plot")

@attr(ANSWER_TEST_TAG)
def test_multi_line_plot():
    ds = fake_random_ds(4)
    fields = [field for field in ds.field_list if field[0] == 'stream']
    field_labels = {f: f[1] for f in fields}
    lines = []
    lines.append(yt.LineBuffer(ds, [0.25, 0, 0], [0.25, 1, 0], 100,
                               label='x = 0.5'))
    lines.append(yt.LineBuffer(ds, [0.5, 0, 0], [0.5, 1, 0], 100,
                               label='x = 0.5'))
    plot = yt.LinePlot.from_lines(ds, fields, lines, field_labels=field_labels)
    plot.annotate_legend(fields[0])
    plot.annotate_legend(fields[1])
    yield compare(ds, plot, test_prefix="answers_multi_line_plot",
                  test_name="answers_multi_line_plot")

def test_line_buffer():
    ds = fake_random_ds(32)
    lb = yt.LineBuffer(ds, (0, 0, 0), (1, 1, 1), 512, label='diag')
    lb['density']
    lb['velocity_x']
    assert_equal(lb['density'].size, 512)
    lb['density'] = 0
    assert_equal(lb['density'], 0)
    assert_equal(set(lb.keys()), set(['density', 'velocity_x']))
    del lb['velocity_x']
    assert_equal(set(lb.keys()), set(['density']))

def test_validate_point():
    ds = fake_random_ds(3)
    with assert_raises(RuntimeError) as ex:
        _validate_point(0, ds, start=True)
    assert_equal(str(ex.exception), 'Input point must be array-like')

    with assert_raises(RuntimeError) as ex:
        _validate_point(ds.arr([[0], [1]], 'code_length'), ds, start=True)
    assert_equal(str(ex.exception), 'Input point must be a 1D array')

    with assert_raises(RuntimeError) as ex:
        _validate_point(ds.arr([0, 1], 'code_length'), ds, start=True)
    assert_equal(str(ex.exception),
                 'Input point must have an element for each dimension')

    ds = fake_random_ds([32, 32, 1])
    _validate_point(ds.arr([0, 1], 'code_length'), ds, start=True)
    _validate_point(ds.arr([0, 1], 'code_length'), ds)
