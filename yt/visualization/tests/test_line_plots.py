import pytest
from numpy.testing import assert_equal

import yt
from yt.testing import fake_random_ds
from yt.visualization.line_plot import _validate_point


def setup_module():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class TestLinePlotSimple:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_random_ds(4)
        fields = [field for field in cls.ds.field_list if field[0] == "stream"]
        field_labels = {f: f[1] for f in fields}
        plot = yt.LinePlot(
            cls.ds, fields, (0, 0, 0), (1, 1, 0), 1000, field_labels=field_labels
        )
        plot.annotate_legend(fields[0])
        plot.annotate_legend(fields[1])
        plot.set_x_unit("cm")
        plot.set_unit(fields[0], "kg/cm**3")
        plot.annotate_title(fields[0], "Density Plot")
        plot.render()
        cls.plot = plot

    @pytest.mark.mpl_image_compare
    def test_lineplot_simple_density(self):
        return self.plot.plots["stream", "density"].figure

    @pytest.mark.mpl_image_compare
    def test_lineplot_simple_velocity_x(self):
        return self.plot.plots["stream", "velocity_x"].figure


class TestLinePlotMulti:
    @classmethod
    def setup_class(cls):
        cls.ds = fake_random_ds(4)
        fields = [field for field in cls.ds.field_list if field[0] == "stream"]
        field_labels = {f: f[1] for f in fields}
        lines = []
        lines.append(
            yt.LineBuffer(cls.ds, [0.25, 0, 0], [0.25, 1, 0], 100, label="x = 0.5")
        )
        lines.append(
            yt.LineBuffer(cls.ds, [0.5, 0, 0], [0.5, 1, 0], 100, label="x = 0.5")
        )
        plot = yt.LinePlot.from_lines(cls.ds, fields, lines, field_labels=field_labels)
        plot.annotate_legend(fields[0])
        plot.annotate_legend(fields[1])
        plot.set_x_unit("cm")
        plot.set_unit(fields[0], "kg/cm**3")
        plot.annotate_title(fields[0], "Density Plot")
        plot.render()
        cls.plot = plot

    @pytest.mark.mpl_image_compare
    def test_lineplot_multi_density(self):
        return self.plot.plots["stream", "density"].figure

    @pytest.mark.mpl_image_compare
    def test_lineplot_multi_velocity_x(self):
        return self.plot.plots["stream", "velocity_x"].figure


def test_line_buffer():
    ds = fake_random_ds(32)
    lb = yt.LineBuffer(ds, (0, 0, 0), (1, 1, 1), 512, label="diag")
    lb[("gas", "density")]
    lb[("gas", "velocity_x")]
    assert_equal(lb[("gas", "density")].size, 512)
    lb[("gas", "density")] = 0
    assert_equal(lb[("gas", "density")], 0)
    assert_equal(set(lb.keys()), {("gas", "density"), ("gas", "velocity_x")})
    del lb[("gas", "velocity_x")]
    assert_equal(set(lb.keys()), {("gas", "density")})


def test_validate_point():
    ds = fake_random_ds(3)
    with pytest.raises(RuntimeError, match="Input point must be array-like"):
        _validate_point(0, ds, start=True)

    with pytest.raises(RuntimeError, match="Input point must be a 1D array"):
        _validate_point(ds.arr([[0], [1]], "code_length"), ds, start=True)

    with pytest.raises(
        RuntimeError, match="Input point must have an element for each dimension"
    ):
        _validate_point(ds.arr([0, 1], "code_length"), ds, start=True)

    ds = fake_random_ds([32, 32, 1])
    _validate_point(ds.arr([0, 1], "code_length"), ds, start=True)
    _validate_point(ds.arr([0, 1], "code_length"), ds)
