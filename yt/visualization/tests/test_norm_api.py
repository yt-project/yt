import matplotlib
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from nose.plugins.attrib import attr
from packaging.version import Version

from yt.testing import (
    ANSWER_TEST_TAG,
    add_noise_fields,
    fake_particle_ds,
    fake_random_ds,
    skip_case,
)
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import (
    LinePlot,
    ParticleProjectionPlot,
    PhasePlot,
    ProfilePlot,
    SlicePlot,
)

ds = fake_random_ds(16)
add_noise_fields(ds)
noise_fields = ["noise0", "noise1", "noise2", "noise3"]
MPL_VERSION = Version(matplotlib.__version__)


@attr(ANSWER_TEST_TAG)
def test_lineplot_set_axis_properties():
    def create_image(filename_prefix):
        p = LinePlot(
            ds,
            ("gas", "density"),
            start_point=[0, 0, 0],
            end_point=[1, 1, 1],
            npoints=32,
        )
        p.set_x_unit("cm")
        p.save(f"{filename_prefix}_xunit")

        p.set_unit(("gas", "density"), "kg/cm**3")
        p.save(f"{filename_prefix}_xunit_zunit")

        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_xunit_zunit_lin")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_lineplot_set_axis_properties"
    test_lineplot_set_axis_properties.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_profileplot_set_axis_properties():
    def create_image(filename_prefix):
        disk = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], (10, "m"), (1, "m"))
        p = ProfilePlot(disk, ("gas", "density"), [("gas", "velocity_x")])
        p.save(f"{filename_prefix}_defaults")

        p.set_unit(("gas", "density"), "kg/cm**3")
        p.save(f"{filename_prefix}_xunit")

        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_xunit_xlin")

        p.set_unit(("gas", "velocity_x"), "mile/hour")
        p.save(f"{filename_prefix}_xunit_xlin_yunit")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_profileplot_set_axis_properties"
    test_profileplot_set_axis_properties.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_sliceplot_custom_norm():
    norms_to_test = [
        (Normalize(), "linear"),
        (LogNorm(), "log"),
    ]
    if MPL_VERSION >= Version("3.4"):
        # don't test this in older versions of MPL because of a deprecation warning
        # that cannot be easily filtered
        norms_to_test.append((SymLogNorm(linthresh=0.01, vmin=-1, vmax=1), "symlog"))

    def create_image(filename_prefix):
        field = ("gas", "density")
        for norm, name in norms_to_test:
            p = SlicePlot(ds, "z", field)
            p.set_norm(field, norm=norm)
            p.save(f"{filename_prefix}_{name}")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_custom_norm"
    test_sliceplot_custom_norm.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_two_slope_norm():
    if MPL_VERSION < Version("3.2"):
        skip_case(reason="TwoSlopeNorm is new in matplotlib 3.2")
    from matplotlib.colors import TwoSlopeNorm

    def create_image(filename_prefix):
        p = SlicePlot(ds, "z", noise_fields)
        p.set_norm("all", norm=TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1))
        # use a diverging colormap for legibility
        p.set_cmap("all", "RdBu")
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_two_slope_norm"
    test_two_slope_norm.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_log():
    def create_image(filename_prefix):
        p = SlicePlot(ds, "z", noise_fields)
        p.set_log("all", False)
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_log"
    test_sliceplot_set_log.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_zlim_and_unit():
    def create_image(filename_prefix):
        field = ("gas", "density")
        p = SlicePlot(ds, "z", field)
        p.set_zlim(field, zmin=0)
        p.set_unit(field, "kg/m**3")
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_zlim_and_unit"
    test_sliceplot_set_zlim_and_unit.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_unit_and_zlim():
    def create_image(filename_prefix):
        field = ("gas", "density")
        p = SlicePlot(ds, "z", field)
        p.set_unit(field, "kg/m**3")
        p.set_zlim(field, zmin=0)
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_unit_and_zlim"
    test_sliceplot_set_unit_and_zlim.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_sliceplot_set_background_color():
    # see https://github.com/yt-project/yt/issues/3854
    def create_image(filename_prefix):
        field = ("gas", "density")
        p = SlicePlot(ds, "z", field, width=1.5)
        p.set_background_color(field, color="C0")
        p.save(f"{filename_prefix}_log")
        p.set_log(("gas", "density"), False)
        p.save(f"{filename_prefix}_lin")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_set_background_color"
    test_sliceplot_set_background_color.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_phaseplot_set_colorbar_properties_implicit():
    def create_image(filename_prefix):
        my_sphere = ds.sphere("c", 1)
        p = PhasePlot(
            my_sphere,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        # using implicit units
        p.set_zlim(("gas", "density"), zmax=10)
        p.save(f"{filename_prefix}_set_zlim_implicit")

        # changing units should affect the colorbar and not the image
        p.set_unit(("gas", "density"), "kg/AU**3")
        p.save(f"{filename_prefix}_set_zlim_set_unit_implicit")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_phaseplot_set_colorbar_properties_implicit"
    test_phaseplot_set_colorbar_properties_implicit.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_phaseplot_set_colorbar_properties_explicit():
    def create_image(filename_prefix):
        my_sphere = ds.sphere("c", 1)
        p = PhasePlot(
            my_sphere,
            ("gas", "noise1"),
            ("gas", "noise3"),
            [("gas", "density")],
            weight_field=None,
        )
        # using explicit units, we expect the colorbar units to stay unchanged
        p.set_zlim(("gas", "density"), zmin=(1e36, "kg/AU**3"))
        p.save(f"{filename_prefix}_set_zlim_explicit")

        # ... until we set them explicitly
        p.set_unit(("gas", "density"), "kg/AU**3")
        p.save(f"{filename_prefix}_set_zlim_set_unit_explicit")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_phaseplot_set_colorbar_properties_explicit"
    test_phaseplot_set_colorbar_properties_explicit.__name__ = test.description
    yield test


@attr(ANSWER_TEST_TAG)
def test_particleprojectionplot_set_colorbar_properties():
    def create_image(filename_prefix):
        field = ("all", "particle_mass")
        p = ParticleProjectionPlot(fake_particle_ds(npart=100), 2, field)
        p.set_buff_size(10)

        p.set_unit(field, "Msun")
        p.save(f"{filename_prefix}_set_unit")

        p.set_zlim(field, zmax=1e-35)
        p.save(f"{filename_prefix}_set_unit_zlim")

        p.set_log(field, False)
        p.save(f"{filename_prefix}_set_unit_zlim_log")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_particleprojectionplot_set_colorbar_properties"
    test_particleprojectionplot_set_colorbar_properties.__name__ = test.description
    yield test
