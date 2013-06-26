"""
Tests for AMRSlice

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Author: Kacper Kowalik <xarthisius.kk@gmail.com>
Affiliation: CA UMK
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Samuel Skillman.  All Rights Reserved.
  Copyright (C) 2013 Kacper Kowalik.  All Rights Reserved.

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
import numpy as np
from nose.tools import raises
from yt.testing import \
    fake_random_pf, assert_equal, assert_array_equal
from yt.utilities.definitions import \
    x_dict, y_dict
from yt.utilities.exceptions import \
    YTNoDataInObjectError

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def teardown_func(fns):
    for fn in fns:
        os.remove(fn)


def test_slice():
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        pf = fake_random_pf(64, nprocs=nprocs)
        dims = pf.domain_dimensions
        xn, yn, zn = pf.domain_dimensions
        xi, yi, zi = pf.domain_left_edge + 1.0 / (pf.domain_dimensions * 2)
        xf, yf, zf = pf.domain_right_edge - 1.0 / (pf.domain_dimensions * 2)
        coords = np.mgrid[xi:xf:xn * 1j, yi:yf:yn * 1j, zi:zf:zn * 1j]
        uc = [np.unique(c) for c in coords]
        slc_pos = 0.5
        # Some simple slice tests with single grids
        for ax, an in enumerate("xyz"):
            xax = x_dict[ax]
            yax = y_dict[ax]
            for wf in ["density", None]:
                fns = []
                slc = pf.h.slice(ax, slc_pos)
                yield assert_equal, slc["ones"].sum(), slc["ones"].size
                yield assert_equal, slc["ones"].min(), 1.0
                yield assert_equal, slc["ones"].max(), 1.0
                yield assert_equal, np.unique(slc["px"]), uc[xax]
                yield assert_equal, np.unique(slc["py"]), uc[yax]
                yield assert_equal, np.unique(slc["pdx"]), 0.5 / dims[xax]
                yield assert_equal, np.unique(slc["pdy"]), 0.5 / dims[yax]
                pw = slc.to_pw()
                fns += pw.save()
                frb = slc.to_frb((1.0, 'unitary'), 64)
                for slc_field in ['ones', 'density']:
                    yield assert_equal, frb[slc_field].info['data_source'], \
                        slc.__str__()
                    yield assert_equal, frb[slc_field].info['axis'], \
                        ax
                    yield assert_equal, frb[slc_field].info['field'], \
                        slc_field
                    yield assert_equal, frb[slc_field].info['units'], \
                        pf.field_info[slc_field].get_units()
                    yield assert_equal, frb[slc_field].info['xlim'], \
                        frb.bounds[:2]
                    yield assert_equal, frb[slc_field].info['ylim'], \
                        frb.bounds[2:]
                    yield assert_equal, frb[slc_field].info['length_to_cm'], \
                        pf['cm']
                    yield assert_equal, frb[slc_field].info['center'], \
                        slc.center
                    yield assert_equal, frb[slc_field].info['coord'], \
                        slc_pos
                teardown_func(fns)
            # wf == None
            yield assert_equal, wf, None


def test_slice_over_edges():
    pf = fake_random_pf(64, nprocs=8, fields=["density"], negative=[False])
    slc = pf.h.slice(0, 0.0)
    slc["density"]
    slc = pf.h.slice(1, 0.5)
    slc["density"]


def test_slice_over_outer_boundary():
    pf = fake_random_pf(64, nprocs=8, fields=["density"], negative=[False])
    slc = pf.h.slice(2, 1.0)
    slc["density"]
    yield assert_equal, slc["density"].size, 0
