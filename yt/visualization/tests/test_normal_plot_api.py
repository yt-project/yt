from itertools import product

import numpy as np
import pytest

from yt.testing import fake_amr_ds
from yt.visualization.plot_window import ProjectionPlot, SlicePlot


@pytest.fixture(scope="module")
def ds():
    return fake_amr_ds(geometry="cartesian")


@pytest.mark.parametrize("plot_cls", (SlicePlot, ProjectionPlot))
def test_normalplot_all_positional_args(ds, plot_cls):
    plot_cls(ds, "z", ("stream", "Density"))


@pytest.mark.parametrize("plot_cls", (SlicePlot, ProjectionPlot))
def test_normalplot_normal_kwarg(ds, plot_cls):
    plot_cls(ds, normal="z", fields=("stream", "Density"))


@pytest.mark.parametrize("plot_cls", (SlicePlot, ProjectionPlot))
def test_error_with_missing_fields_and_normal(ds, plot_cls):
    with pytest.raises(
        TypeError,
        match="missing 2 required positional arguments: 'normal' and 'fields'",
    ):
        plot_cls(ds)


@pytest.mark.parametrize("plot_cls", (SlicePlot, ProjectionPlot))
def test_error_with_missing_fields_with_normal_kwarg(ds, plot_cls):
    with pytest.raises(
        TypeError, match=r"missing (1 )?required positional argument: 'fields'$"
    ):
        plot_cls(ds, normal="z")


@pytest.mark.parametrize("plot_cls", (SlicePlot, ProjectionPlot))
def test_error_with_missing_fields_with_positional(ds, plot_cls):
    with pytest.raises(
        TypeError, match=r"missing (1 )?required positional argument: 'fields'$"
    ):
        plot_cls(ds, "z")


@pytest.mark.parametrize(
    "plot_cls, normal",
    product([SlicePlot, ProjectionPlot], [(0, 0, 1), [0, 0, 1], np.array((0, 0, 1))]),
)
def test_normalplot_normal_array(ds, plot_cls, normal):
    # see regression https://github.com/yt-project/yt/issues/3736
    plot_cls(ds, normal, fields=("stream", "Density"))
