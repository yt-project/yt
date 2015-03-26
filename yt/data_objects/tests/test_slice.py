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
    fake_random_ds, assert_equal
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
    grid_eps = np.finfo(np.float64).eps
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        ds = fake_random_ds(64, nprocs=nprocs)
        dims = ds.domain_dimensions
        xn, yn, zn = ds.domain_dimensions
        dx = ds.arr(1.0 / (ds.domain_dimensions * 2), 'code_length')
        xi, yi, zi = ds.domain_left_edge + dx
        xf, yf, zf = ds.domain_right_edge - dx
        coords = np.mgrid[xi:xf:xn * 1j, yi:yf:yn * 1j, zi:zf:zn * 1j]
        uc = [np.unique(c) for c in coords]
        slc_pos = 0.5
        # Some simple slice tests with single grids
        for ax, an in enumerate("xyz"):
            xax = ds.coordinates.x_axis[ax]
            yax = ds.coordinates.y_axis[ax]
            for wf in ["density", None]:
                slc = ds.slice(ax, slc_pos)
                shifted_slc = ds.slice(ax, slc_pos + grid_eps)
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
                for width in [(1.0, 'unitary'), 1.0, ds.quan(0.5, 'code_length')]:
                    frb = slc.to_frb((1.0, 'unitary'), 64)
                    shifted_frb = shifted_slc.to_frb((1.0, 'unitary'), 64)
                    for slc_field in ['ones', 'density']:
                        fi = ds._get_field_info(slc_field)
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
                        yield assert_equal, frb[slc_field], \
                            shifted_frb[slc_field]
            yield assert_equal, wf, None
    teardown_func(fns)


def test_slice_over_edges():
    ds = fake_random_ds(64, nprocs=8, fields=["density"], negative=[False])
    slc = ds.slice(0, 0.0)
    slc["density"]
    slc = ds.slice(1, 0.5)
    slc["density"]


def test_slice_over_outer_boundary():
    ds = fake_random_ds(64, nprocs=8, fields=["density"], negative=[False])
    slc = ds.slice(2, 1.0)
    slc["density"]
    yield assert_equal, slc["density"].size, 0
