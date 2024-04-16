import numpy as np
import pytest

from yt import SlicePlot, load_uniform_grid
from yt.testing import fake_amr_ds, requires_module


@requires_module("cartopy")
@pytest.mark.parametrize("geometry", ["geographic", "internal_geographic"])
def test_quiver_callback_geographic(geometry):
    flds = ("density", "velocity_ew", "velocity_ns")
    units = ("g/cm**3", "m/s", "m/s")

    ds = fake_amr_ds(fields=flds, units=units, geometry=geometry)

    for ax in ds.coordinates.axis_order:
        slc = SlicePlot(ds, ax, "density", buff_size=(50, 50))
        if ax == ds.coordinates.radial_axis:
            # avoid the exact transform bounds
            slc.set_width((359.99, 179.99))
        slc.annotate_quiver(("stream", "velocity_ew"), ("stream", "velocity_ns"))
        slc.render()


@pytest.fixture()
def ds_geo_uni_grid():
    yc = 0.0
    xc = 0.0

    def _vel_calculator(grid, ax):
        y_lat = grid.fcoords[:, 1].d
        x_lon = grid.fcoords[:, 2].d
        x_lon[x_lon > 180] = x_lon[x_lon > 180] - 360.0
        dist = np.sqrt((y_lat - yc) ** 2 + (xc - x_lon) ** 2)
        if ax == 1:
            sgn = np.sign(y_lat - yc)
        elif ax == 2:
            sgn = np.sign(x_lon - xc)

        vel = np.exp(-((dist / 45) ** 2)) * sgn
        return vel.reshape(grid.shape)

    def _calculate_u(grid, field):
        return _vel_calculator(grid, 2)

    def _calculate_v(grid, field):
        return _vel_calculator(grid, 1)

    data = {"u_vel": _calculate_u, "v_vel": _calculate_v}

    bbox = [[10.0, 1000], [-90, 90], [0, 360]]
    bbox = np.array(bbox)

    ds = load_uniform_grid(
        data,
        (16, 16, 32),
        bbox=bbox,
        geometry="geographic",
        axis_order=("altitude", "latitude", "longitude"),
    )

    return ds


@requires_module("cartopy")
@pytest.mark.mpl_image_compare
def test_geoquiver_answer(ds_geo_uni_grid):
    slc = SlicePlot(ds_geo_uni_grid, "altitude", "u_vel")
    slc.set_width((359.99, 179.99))
    slc.set_log("u_vel", False)
    slc.annotate_quiver("u_vel", "v_vel", scale=50)
    slc.render()
    return slc.plots["u_vel"].figure
