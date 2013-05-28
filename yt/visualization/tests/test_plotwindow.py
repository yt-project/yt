"""
Testsuite for PlotWindow class

Author: Nathan Goldbaum <goldbaum@ucolick.org>
Affiliation: UCSC Astronomy
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Nathan Goldbaum.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import os
import tempfile
import shutil
from yt.testing import \
    fake_random_pf, assert_equal, assert_rel_equal
from yt.utilities.answer_testing.framework import \
    requires_pf, data_dir_load, PlotWindowAttributeTest
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    with open(fname, 'rb') as fimg:
        data = fimg.read()
    data = str(data)
    image_type = ''

    # see http://www.w3.org/TR/PNG/#5PNG-file-signature
    if data.startswith('\211PNG\r\n\032\n'):
        image_type = '.png'
    # see http://www.mathguide.de/info/tools/media-types/image/jpeg
    elif data.startswith('\377\330'):
        image_type = '.jpeg'
    elif data.startswith('%!PS-Adobe'):
        if 'EPSF' in data[:data.index('\n')]:
            image_type = '.eps'
        else:
            image_type = '.ps'
    elif data.startswith('%PDF'):
        image_type = '.pdf'

    return image_type == os.path.splitext(fname)[1]

attr_args ={ "pan"             : [( ((0.1, 0.1),), {} )],
             "pan_rel"         : [( ((0.1, 0.1),), {} )],
             "set_axes_unit"   : [( ("kpc",), {} ),
                                  ( ("Mpc",), {} ),
                                  ( (("kpc", "kpc"),), {} ),
                                  ( (("kpc", "Mpc"),), {} )],
             "set_buff_size"   : [( (1600,), {} ),
                                  ( ((600, 800),), {} )],
             "set_center"      : [( ((0.4, 0.3),), {} )],
             "set_cmap"        : [( ('Density', 'RdBu'), {} ),
                                  ( ('Density', 'kamae'), {} )],
             "set_font"        : [( ({'family':'sans-serif', 'style':'italic',
                                      'weight':'bold', 'size':24},), {} )],
             "set_log"         : [( ('Density', False), {} )],
             "set_window_size" : [( (7.0,), {} )],
             "set_zlim" : [( ('Density', 1e-25, 1e-23), {} ),
                           ( ('Density', 1e-25, None), {'dynamic_range' : 4} )],
             "zoom" : [( (10,), {} )] }

m7 = "DD0010/moving7_0010"
wt = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"
@requires_pf(m7)
@requires_pf(wt)
def test_attributes():
    """Test plot member functions that aren't callbacks"""
    plot_field = 'Density'
    decimals = 3

    pf = data_dir_load(m7)
    for ax in 'xyz':
        for attr_name in attr_args.keys():
            for args in attr_args[attr_name]:
                yield PlotWindowAttributeTest(pf, plot_field, ax, attr_name,
                                              args, decimals)
    pf = data_dir_load(wt)
    ax = 'z'
    for attr_name in attr_args.keys():
        for args in attr_args[attr_name]:
            yield PlotWindowAttributeTest(pf, plot_field, ax, attr_name,
                                          args, decimals)

def test_setwidth():
    pf = fake_random_pf(64)

    slc = SlicePlot(pf, 0, 'Density')

    yield assert_equal, [slc.xlim, slc.ylim, slc.width], \
        [(0.0, 1.0), (0.0, 1.0), (1.0, 1.0)]

    slc.set_width((0.5,0.8))

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(0.25, 0.75), (0.1, 0.9), (0.5, 0.8)], 15

    slc.set_width(15,'kpc')

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(-7.5/pf['kpc'], 7.5/pf['kpc']),
         (-7.5/pf['kpc'], 7.5/pf['kpc']),
         (15/pf['kpc'], 15/pf['kpc'])], 15

    slc.set_width((15,'kpc'))

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(-7.5/pf['kpc'], 7.5/pf['kpc']),
         (-7.5/pf['kpc'], 7.5/pf['kpc']),
         (15/pf['kpc'], 15/pf['kpc'])], 15

    slc.set_width(((15,'kpc'),(10,'kpc')))

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(-7.5/pf['kpc'], 7.5/pf['kpc']),
         (-5/pf['kpc'], 5/pf['kpc']),
         (15/pf['kpc'], 10/pf['kpc'])], 15

    slc.set_width(((15,'kpc'),(10000,'pc')))

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(-7.5/pf['kpc'], 7.5/pf['kpc']),
         (-5/pf['kpc'], 5/pf['kpc']),
         (15/pf['kpc'], 10/pf['kpc'])], 15

    slc.set_width((15,'kpc'),(10000,'pc'))

    yield assert_rel_equal, [slc.xlim, slc.ylim, slc.width], \
        [(-7.5/pf['kpc'], 7.5/pf['kpc']),
         (-5/pf['kpc'], 5/pf['kpc']),
         (15/pf['kpc'], 10/pf['kpc'])], 15

def test_save():
    """Test plot window creation and saving to disk."""
    # Perform I/O in safe place instead of yt main dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    normal = [1, 1, 1]

    test_pf = fake_random_pf(64)
    test_flnms = [None, 'test.png', 'test.eps',
                  'test.ps', 'test.pdf']

    ds_region = test_pf.h.region([0.5]*3,[0.4]*3,[0.6]*3)

    for dim in [0, 1, 2]:
        obj = SlicePlot(test_pf, dim, 'Density')
        for fname in test_flnms:
            yield assert_equal, assert_fname(obj.save(fname)[0]), True

    for dim in [0, 1, 2]:
        obj = ProjectionPlot(test_pf, dim, 'Density')
        for fname in test_flnms:
            yield assert_equal, assert_fname(obj.save(fname)[0]), True
        # Test ProjectionPlot's data_source keyword
        obj = ProjectionPlot(test_pf, dim, 'Density',
                             data_source=ds_region)
        obj.save()

    obj = OffAxisSlicePlot(test_pf, normal, 'Density')
    for fname in test_flnms:
        yield assert_equal, assert_fname(obj.save(fname)[0]), True

    obj = OffAxisProjectionPlot(test_pf, normal, 'Density')
    for fname in test_flnms:
        yield assert_equal, assert_fname(obj.save(fname)[0]), True

    os.chdir(curdir)
    # clean up
    shutil.rmtree(tmpdir)
