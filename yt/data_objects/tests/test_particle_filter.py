import os
import shutil
import tempfile

import numpy as np
from nose.tools import assert_raises

from yt.data_objects.particle_filters import add_particle_filter, particle_filter
from yt.testing import assert_equal, fake_random_ds, fake_sph_grid_ds
from yt.utilities.exceptions import YTIllDefinedFilter, YTIllDefinedParticleFilter
from yt.visualization.plot_window import ProjectionPlot


def test_add_particle_filter():
    """Test particle filters created via add_particle_filter

    This accesses a deposition field using the particle filter, which was a
    problem in previous versions on this dataset because there are chunks with
    no stars in them.

    """

    def stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] > 0.5

    add_particle_filter(
        "stars1", function=stars, filtered_type="all", requires=["particle_mass"]
    )
    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("stars1")
    assert ("deposit", "stars1_cic") in ds.derived_field_list

    # Test without requires field
    add_particle_filter("stars2", function=stars)
    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("stars2")
    assert ("deposit", "stars2_cic") in ds.derived_field_list

    # Test adding filter with fields not defined on the ds
    with assert_raises(YTIllDefinedParticleFilter) as ex:
        add_particle_filter(
            "bad_stars", function=stars, filtered_type="all", requires=["wrong_field"]
        )
        ds.add_particle_filter("bad_stars")
    actual = str(ex.exception)
    desired = (
        "\nThe fields\n\t('all', 'wrong_field'),\nrequired by the"
        ' "bad_stars" particle filter, are not defined for this dataset.'
    )
    assert_equal(actual, desired)


def test_add_particle_filter_overriding():
    """Test the add_particle_filter overriding"""
    from yt.data_objects.particle_filters import filter_registry
    from yt.funcs import mylog

    def star_0(pfilter, data):
        pass

    def star_1(pfilter, data):
        pass

    # Use a closure to store whether the warning was called
    def closure(status):
        def warning_patch(*args, **kwargs):
            status[0] = True

        def was_called():
            return status[0]

        return warning_patch, was_called

    ## Test 1: we add a dummy particle filter
    add_particle_filter(
        "dummy", function=star_0, filtered_type="all", requires=["creation_time"]
    )
    assert "dummy" in filter_registry
    assert_equal(filter_registry["dummy"].function, star_0)

    ## Test 2: we add another dummy particle filter.
    ##         a warning is expected. We use the above closure to
    ##         check that.
    # Store the original warning function
    warning = mylog.warning
    monkey_warning, monkey_patch_was_called = closure([False])
    mylog.warning = monkey_warning
    add_particle_filter(
        "dummy", function=star_1, filtered_type="all", requires=["creation_time"]
    )
    assert_equal(filter_registry["dummy"].function, star_1)
    assert_equal(monkey_patch_was_called(), True)

    # Restore the original warning function
    mylog.warning = warning


def test_particle_filter_decorator():
    """Test the particle_filter decorator"""

    @particle_filter(filtered_type="all", requires=["particle_mass"])
    def heavy_stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] > 0.5

    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("heavy_stars")
    assert "heavy_stars" in ds.particle_types
    assert ("deposit", "heavy_stars_cic") in ds.derived_field_list

    # Test name of particle filter
    @particle_filter(name="my_stars", filtered_type="all", requires=["particle_mass"])
    def custom_stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] == 0.5

    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("my_stars")
    assert "my_stars" in ds.particle_types
    assert ("deposit", "my_stars_cic") in ds.derived_field_list


def test_particle_filter_exceptions():
    @particle_filter(filtered_type="all", requires=["particle_mass"])
    def filter1(pfilter, data):
        return data

    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("filter1")

    ad = ds.all_data()
    with assert_raises(YTIllDefinedFilter):
        ad["filter1", "particle_mass"].shape[0]

    @particle_filter(filtered_type="all", requires=["particle_mass"])
    def filter2(pfilter, data):
        filter_field = ("io", "particle_mass")
        return data[filter_field] > 0.5

    ds.add_particle_filter("filter2")
    ad = ds.all_data()
    ad["filter2", "particle_mass"].min()


def test_particle_filter_dependency():
    """
    Test dataset add_particle_filter which should automatically add
    the dependency of the filter.
    """

    @particle_filter(filtered_type="all", requires=["particle_mass"])
    def h_stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] > 0.5

    @particle_filter(filtered_type="h_stars", requires=["particle_mass"])
    def hh_stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] > 0.9

    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("hh_stars")
    assert "hh_stars" in ds.particle_types
    assert "h_stars" in ds.particle_types
    assert ("deposit", "hh_stars_cic") in ds.derived_field_list
    assert ("deposit", "h_stars_cic") in ds.derived_field_list


def test_covering_grid_particle_filter():
    @particle_filter(filtered_type="all", requires=["particle_mass"])
    def heavy_stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "particle_mass")
        return data[filter_field] > 0.5

    ds = fake_random_ds(16, nprocs=8, particles=16)
    ds.add_particle_filter("heavy_stars")

    for grid in ds.index.grids:
        cg = ds.covering_grid(grid.Level, grid.LeftEdge, grid.ActiveDimensions)

        assert_equal(
            cg["heavy_stars", "particle_mass"].shape[0],
            grid["heavy_stars", "particle_mass"].shape[0],
        )
        assert_equal(
            cg["heavy_stars", "particle_mass"].shape[0],
            grid["heavy_stars", "particle_mass"].shape[0],
        )


def test_sph_particle_filter_plotting():
    ds = fake_sph_grid_ds()

    @particle_filter("central_gas", requires=["particle_position"], filtered_type="io")
    def _filter(pfilter, data):
        coords = np.abs(data[pfilter.filtered_type, "particle_position"])
        return (coords[:, 0] < 1.6) & (coords[:, 1] < 1.6) & (coords[:, 2] < 1.6)

    ds.add_particle_filter("central_gas")

    plot = ProjectionPlot(ds, "z", ("central_gas", "density"))
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    plot.save()

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
