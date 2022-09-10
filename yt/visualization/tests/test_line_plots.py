import pytest

import yt
from yt.testing import assert_equal, fake_random_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.line_plot import _validate_point


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def compare(ds, plot, test_prefix, test_name, decimals=12):
    def image_from_plot(filename_prefix):
        return plot.save(filename_prefix)

    image_from_plot.__name__ = f"line_{test_prefix}"
    test = GenericImageTest(ds, image_from_plot, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test


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
        cls.plot = plot

    @pytest.mark.mpl_image_compare(filename="line_plot_simple_density.png")
    def test_line_plot_simple_density(self):
        return self.plot.plots["stream", "density"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_simple_velocity_x.png")
    def test_line_plot_simple_velocity_x(self):
        return self.plot.plots["stream", "velocity_x"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_simple_velocity_y.png")
    def test_line_plot_simple_velocity_y(self):
        return self.plot.plots["stream", "velocity_y"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_simple_velocity_z.png")
    def test_line_plot_simple_velocity_z(self):
        return self.plot.plots["stream", "velocity_z"].figure


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
        cls.plot = plot

    @pytest.mark.mpl_image_compare(filename="line_plot_multi_density.png")
    def test_line_plot_multi_density(self):
        return self.plot.plots["stream", "density"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_multi_velocity_x.png")
    def test_line_plot_multi_velocity_x(self):
        return self.plot.plots["stream", "velocity_x"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_multi_velocity_y.png")
    def test_line_plot_multi_velocity_y(self):
        return self.plot.plots["stream", "velocity_y"].figure

    @pytest.mark.mpl_image_compare(filename="line_plot_multi_velocity_z.png")
    def test_line_plot_multi_velocity_z(self):
        return self.plot.plots["stream", "velocity_z"].figure


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
