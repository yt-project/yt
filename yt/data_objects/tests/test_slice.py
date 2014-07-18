"""
Tests for AMRSlice


"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import os
import numpy as np
import tempfile
from yt.testing import \
    fake_random_pf, assert_equal
from yt.units.unit_object import Unit


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def teardown_func(fns):
    for fn in fns:
        try:
            os.remove(fn)
        except OSError:
            pass


def test_slice():
    fns = []
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        pf = fake_random_pf(64, nprocs=nprocs)
        dims = pf.domain_dimensions
        xn, yn, zn = pf.domain_dimensions
        dx = pf.arr(1.0 / (pf.domain_dimensions * 2), 'code_length')
        xi, yi, zi = pf.domain_left_edge + dx
        xf, yf, zf = pf.domain_right_edge - dx
        coords = np.mgrid[xi:xf:xn * 1j, yi:yf:yn * 1j, zi:zf:zn * 1j]
        uc = [np.unique(c) for c in coords]
        slc_pos = 0.5
        # Some simple slice tests with single grids
        for ax, an in enumerate("xyz"):
            xax = pf.coordinates.x_axis[ax]
            yax = pf.coordinates.y_axis[ax]
            for wf in ["density", None]:
                slc = pf.slice(ax, slc_pos)
                yield assert_equal, slc["ones"].sum(), slc["ones"].size
                yield assert_equal, slc["ones"].min(), 1.0
                yield assert_equal, slc["ones"].max(), 1.0
                yield assert_equal, np.unique(slc["px"]), uc[xax]
                yield assert_equal, np.unique(slc["py"]), uc[yax]
                yield assert_equal, np.unique(slc["pdx"]), 0.5 / dims[xax]
                yield assert_equal, np.unique(slc["pdy"]), 0.5 / dims[yax]
                pw = slc.to_pw(fields='density')
                for p in pw.plots.values():
                    tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
                    os.close(tmpfd)
                    p.save(name=tmpname)
                    fns.append(tmpname)
                frb = slc.to_frb((1.0, 'unitary'), 64)
                for slc_field in ['ones', 'density']:
                    fi = pf._get_field_info(slc_field)
                    yield assert_equal, frb[slc_field].info['data_source'], \
                        slc.__str__()
                    yield assert_equal, frb[slc_field].info['axis'], \
                        ax
                    yield assert_equal, frb[slc_field].info['field'], \
                        slc_field
                    yield assert_equal, frb[slc_field].units, \
                        Unit(fi.units)
                    yield assert_equal, frb[slc_field].info['xlim'], \
                        frb.bounds[:2]
                    yield assert_equal, frb[slc_field].info['ylim'], \
                        frb.bounds[2:]
                    yield assert_equal, frb[slc_field].info['center'], \
                        slc.center
                    yield assert_equal, frb[slc_field].info['coord'], \
                        slc_pos
            # wf == None
            yield assert_equal, wf, None
    teardown_func(fns)


def test_slice_over_edges():
    pf = fake_random_pf(64, nprocs=8, fields=["density"], negative=[False])
    slc = pf.slice(0, 0.0)
    slc["density"]
    slc = pf.slice(1, 0.5)
    slc["density"]


def test_slice_over_outer_boundary():
    pf = fake_random_pf(64, nprocs=8, fields=["density"], negative=[False])
    slc = pf.slice(2, 1.0)
    slc["density"]
    yield assert_equal, slc["density"].size, 0
