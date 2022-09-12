import os
import tempfile
from unittest import mock

import numpy as np

from yt.testing import assert_equal, assert_rel_equal, fake_amr_ds, fake_random_ds
from yt.units.unit_object import Unit

LENGTH_UNIT = 2.0


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def teardown_func(fns):
    for fn in fns:
        try:
            os.remove(fn)
        except OSError:
            pass


@mock.patch("yt.visualization._mpl_imports.FigureCanvasAgg.print_figure")
def test_projection(pf):
    fns = []
    for nprocs in [8, 1]:
        # We want to test both 1 proc and 8 procs, to make sure that
        # parallelism isn't broken
        fields = ("density", "temperature", "velocity_x", "velocity_y", "velocity_z")
        units = ("g/cm**3", "K", "cm/s", "cm/s", "cm/s")
        ds = fake_random_ds(
            64, fields=fields, units=units, nprocs=nprocs, length_unit=LENGTH_UNIT
        )
        dims = ds.domain_dimensions
        xn, yn, zn = ds.domain_dimensions
        xi, yi, zi = ds.domain_left_edge.to_ndarray() + 1.0 / (ds.domain_dimensions * 2)
        xf, yf, zf = ds.domain_right_edge.to_ndarray() - 1.0 / (
            ds.domain_dimensions * 2
        )
        dd = ds.all_data()
        coords = np.mgrid[xi : xf : xn * 1j, yi : yf : yn * 1j, zi : zf : zn * 1j]
        uc = [np.unique(c) for c in coords]
        # test if projections inherit the field parameters of their data sources
        dd.set_field_parameter("bulk_velocity", np.array([0, 1, 2]))
        proj = ds.proj(("gas", "density"), 0, data_source=dd)
        assert_equal(
            dd.field_parameters["bulk_velocity"], proj.field_parameters["bulk_velocity"]
        )

        # Some simple projection tests with single grids
        for ax, an in enumerate("xyz"):
            xax = ds.coordinates.x_axis[ax]
            yax = ds.coordinates.y_axis[ax]
            for wf in [("gas", "density"), None]:
                proj = ds.proj(
                    [("index", "ones"), ("gas", "density")], ax, weight_field=wf
                )
                if wf is None:
                    assert_equal(
                        proj[("index", "ones")].sum(),
                        LENGTH_UNIT * proj[("index", "ones")].size,
                    )
                    assert_equal(proj[("index", "ones")].min(), LENGTH_UNIT)
                    assert_equal(proj[("index", "ones")].max(), LENGTH_UNIT)
                else:
                    assert_equal(
                        proj[("index", "ones")].sum(), proj[("index", "ones")].size
                    )
                    assert_equal(proj[("index", "ones")].min(), 1.0)
                    assert_equal(proj[("index", "ones")].max(), 1.0)
                assert_equal(np.unique(proj["px"]), uc[xax])
                assert_equal(np.unique(proj["py"]), uc[yax])
                assert_equal(np.unique(proj["pdx"]), 1.0 / (dims[xax] * 2.0))
                assert_equal(np.unique(proj["pdy"]), 1.0 / (dims[yax] * 2.0))
                plots = [proj.to_pw(fields=("gas", "density")), proj.to_pw()]
                for pw in plots:
                    for p in pw.plots.values():
                        tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
                        os.close(tmpfd)
                        p.save(name=tmpname)
                        fns.append(tmpname)
                frb = proj.to_frb((1.0, "unitary"), 64)
                for proj_field in [
                    ("index", "ones"),
                    ("gas", "density"),
                    ("gas", "temperature"),
                ]:
                    fi = ds._get_field_info(proj_field)
                    assert_equal(frb[proj_field].info["data_source"], proj.__str__())
                    assert_equal(frb[proj_field].info["axis"], ax)
                    assert_equal(frb[proj_field].info["field"], str(proj_field))
                    field_unit = Unit(fi.units)
                    if wf is not None:
                        assert_equal(
                            frb[proj_field].units,
                            Unit(field_unit, registry=ds.unit_registry),
                        )
                    else:
                        if frb[proj_field].units.is_code_unit:
                            proj_unit = "code_length"
                        else:
                            proj_unit = "cm"
                        if field_unit != "" and field_unit != Unit():
                            proj_unit = f"({field_unit}) * {proj_unit}"
                        assert_equal(
                            frb[proj_field].units,
                            Unit(proj_unit, registry=ds.unit_registry),
                        )
                    assert_equal(frb[proj_field].info["xlim"], frb.bounds[:2])
                    assert_equal(frb[proj_field].info["ylim"], frb.bounds[2:])
                    assert_equal(frb[proj_field].info["center"], proj.center)
                    if wf is None:
                        assert_equal(frb[proj_field].info["weight_field"], wf)
                    else:
                        assert_equal(
                            frb[proj_field].info["weight_field"],
                            proj.data_source._determine_fields(wf)[0],
                        )
            # wf == None
            assert_equal(wf, None)
            v1 = proj[("gas", "density")].sum()
            v2 = (dd[("gas", "density")] * dd[("index", f"d{an}")]).sum()
            assert_rel_equal(v1, v2.in_units(v1.units), 10)

        # Test moment projections
        def make_vsq_field(aname):
            def _vsquared(field, data):
                return data["gas", f"velocity_{aname}"] ** 2

            return _vsquared

        for ax, an in enumerate("xyz"):
            ds.add_field(
                ("gas", f"velocity_{an}_squared"),
                make_vsq_field(an),
                sampling_type="local",
                units="cm**2/s**2",
            )
            proj1 = ds.proj(
                [("gas", f"velocity_{an}"), ("gas", f"velocity_{an}_squared")],
                ax,
                weight_field=("gas", "density"),
                moment=1,
            )
            proj2 = ds.proj(
                ("gas", f"velocity_{an}"), ax, weight_field=("gas", "density"), moment=2
            )
            assert_rel_equal(
                np.sqrt(
                    proj1["gas", f"velocity_{an}_squared"]
                    - proj1["gas", f"velocity_{an}"] ** 2
                ),
                proj2["gas", f"velocity_{an}"],
                10,
            )
    teardown_func(fns)


def test_max_level():
    ds = fake_amr_ds(fields=[("gas", "density")], units=["mp/cm**3"])
    proj = ds.proj(("gas", "density"), 2, method="max", max_level=2)
    assert proj[("index", "grid_level")].max() == 2

    proj = ds.proj(("gas", "density"), 2, method="max")
    assert proj[("index", "grid_level")].max() == ds.index.max_level


def test_min_level():
    ds = fake_amr_ds(fields=[("gas", "density")], units=["mp/cm**3"])
    proj = ds.proj(("gas", "density"), 2, method="min")
    assert proj[("index", "grid_level")].min() == 0

    proj = ds.proj(("gas", "density"), 2, method="max")
    assert proj[("index", "grid_level")].min() == ds.index.min_level
