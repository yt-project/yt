import numpy as np
from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_rel_equal, \
    fake_amr_ds
from yt.units.unit_object import Unit
import os
import tempfile

LENGTH_UNIT = 2.0


def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def teardown_func(fns):
    for fn in fns:
        try:
            os.remove(fn)
        except OSError:
            pass

def test_projection():
    fns = []
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        fields = ("density", "temperature", "velocity_x", "velocity_y", "velocity_z")
        units = ('g/cm**3', 'K', 'cm/s', 'cm/s', 'cm/s')
        ds = fake_random_ds(64, fields=fields, units=units, nprocs=nprocs,
                            length_unit=LENGTH_UNIT)
        dims = ds.domain_dimensions
        xn, yn, zn = ds.domain_dimensions
        xi, yi, zi = ds.domain_left_edge.to_ndarray() + \
            1.0 / (ds.domain_dimensions * 2)
        xf, yf, zf = ds.domain_right_edge.to_ndarray() - \
            1.0 / (ds.domain_dimensions * 2)
        dd = ds.all_data()
        coords = np.mgrid[xi:xf:xn*1j, yi:yf:yn*1j, zi:zf:zn*1j]
        uc = [np.unique(c) for c in coords]
        # test if projections inherit the field parameters of their data sources
        dd.set_field_parameter("bulk_velocity", np.array([0,1,2]))
        proj = ds.proj("density", 0, data_source=dd)
        yield assert_equal, dd.field_parameters["bulk_velocity"], \
          proj.field_parameters["bulk_velocity"]

        # Some simple projection tests with single grids
        for ax, an in enumerate("xyz"):
            xax = ds.coordinates.x_axis[ax]
            yax = ds.coordinates.y_axis[ax]
            for wf in ['density', ("gas", "density"), None]:
                proj = ds.proj(["ones", "density"], ax, weight_field=wf)
                if wf is None:
                    yield assert_equal, proj["ones"].sum(), LENGTH_UNIT*proj["ones"].size
                    yield assert_equal, proj["ones"].min(), LENGTH_UNIT
                    yield assert_equal, proj["ones"].max(), LENGTH_UNIT
                else:
                    yield assert_equal, proj["ones"].sum(), proj["ones"].size
                    yield assert_equal, proj["ones"].min(), 1.0
                    yield assert_equal, proj["ones"].max(), 1.0
                yield assert_equal, np.unique(proj["px"]), uc[xax]
                yield assert_equal, np.unique(proj["py"]), uc[yax]
                yield assert_equal, np.unique(proj["pdx"]), 1.0/(dims[xax]*2.0)
                yield assert_equal, np.unique(proj["pdy"]), 1.0/(dims[yax]*2.0)
                plots = [proj.to_pw(fields='density'), proj.to_pw()]
                for pw in plots:
                    for p in pw.plots.values():
                        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
                        os.close(tmpfd)
                        p.save(name=tmpname)
                        fns.append(tmpname)
                frb = proj.to_frb((1.0, 'unitary'), 64)
                for proj_field in ['ones', 'density', 'temperature']:
                    fi = ds._get_field_info(proj_field)
                    yield assert_equal, frb[proj_field].info['data_source'], \
                        proj.__str__()
                    yield assert_equal, frb[proj_field].info['axis'], \
                        ax
                    yield assert_equal, frb[proj_field].info['field'], \
                        proj_field
                    field_unit = Unit(fi.units)
                    if wf is not None:
                        yield assert_equal, frb[proj_field].units, \
                            Unit(field_unit, registry=ds.unit_registry)
                    else:
                        if frb[proj_field].units.is_code_unit:
                            proj_unit = "code_length"
                        else:
                            proj_unit = "cm"
                        if field_unit != '' and field_unit != Unit():
                            proj_unit = \
                                "({0}) * {1}".format(field_unit, proj_unit)
                        yield assert_equal, frb[proj_field].units, \
                            Unit(proj_unit, registry=ds.unit_registry)
                    yield assert_equal, frb[proj_field].info['xlim'], \
                        frb.bounds[:2]
                    yield assert_equal, frb[proj_field].info['ylim'], \
                        frb.bounds[2:]
                    yield assert_equal, frb[proj_field].info['center'], \
                        proj.center
                    if wf is None:
                        yield assert_equal, \
                            frb[proj_field].info['weight_field'], wf
                    else:
                        yield assert_equal, \
                            frb[proj_field].info['weight_field'], \
                            proj.data_source._determine_fields(wf)[0]
            # wf == None
            yield assert_equal, wf, None
            v1 = proj["density"].sum()
            v2 = (LENGTH_UNIT * dd["density"] * dd["d%s" % an]).sum()
            yield assert_rel_equal, v1, v2, 10
    teardown_func(fns)


def test_max_level():
    ds = fake_amr_ds()
    proj = ds.proj('Density', 2, method='mip', max_level=2)
    assert proj['grid_level'].max() == 2

    proj = ds.proj('Density', 2, method='mip')
    assert proj['grid_level'].max() == ds.index.max_level
