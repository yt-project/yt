from __future__ import print_function
import yt
from yt.testing import \
    assert_equal, \
    requires_file
from yt.data_objects.particle_filters import \
    add_particle_filter, particle_filter
import numpy as np


# Dataset required for this test
iso_galaxy = 'IsolatedGalaxy/galaxy0030/galaxy0030'


@requires_file(iso_galaxy)
def test_add_particle_filter():
    """Test particle filters created via add_particle_filter

    This accesses a deposition field using the particle filter, which was a
    problem in previous versions on this dataset because there are chunks with
    no stars in them.

    """

    def stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "creation_time")
        return (data.ds.current_time - data[filter_field]) > 0

    add_particle_filter("stars", function=stars, filtered_type='all',
                        requires=["creation_time"])
    ds = yt.load(iso_galaxy)
    ds.add_particle_filter('stars')
    ad = ds.all_data()
    ad['deposit', 'stars_cic']
    assert True

@requires_file(iso_galaxy)
def test_particle_filter():
    """Test the particle_filter decorator"""

    @particle_filter(filtered_type='all', requires=['creation_time'])
    def stars(pfilter, data):
        filter_field = (pfilter.filtered_type, "creation_time")
        return (data.ds.current_time - data[filter_field]) > 0

    ds = yt.load(iso_galaxy)
    ds.add_particle_filter('stars')
    ad = ds.all_data()
    ad['deposit', 'stars_cic']
    assert True

@requires_file(iso_galaxy)
def test_particle_filter_dependency():
    """
    Test dataset add_particle_filter which should automatically add
    the dependency of the filter.
    """

    @particle_filter(filtered_type='all', requires=['particle_type'])
    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 2
        return filter

    @particle_filter(filtered_type='stars', requires=['creation_time'])
    def young_stars(pfilter, data):
        age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
        filter = np.logical_and(age.in_units('Myr') <= 5, age >= 0)
        return filter

    ds = yt.load(iso_galaxy)
    ds.add_particle_filter('young_stars')
    assert 'young_stars' in ds.particle_types
    assert 'stars' in ds.particle_types
    print(ds.derived_field_list)
    assert ('deposit', 'young_stars_cic') in ds.derived_field_list
    assert ('deposit', 'stars_cic') in ds.derived_field_list

@requires_file(iso_galaxy)
def test_covering_grid_particle_filter():
    @particle_filter(requires=["particle_type"], filtered_type='all')
    def stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 2
        return filter

    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')

    ds.add_particle_filter('stars')

    for grid in ds.index.grids[20:31]:
        cg = ds.covering_grid(grid.Level, grid.LeftEdge, grid.ActiveDimensions)

        assert_equal(cg['stars', 'particle_ones'].shape[0],
                     grid['stars', 'particle_ones'].shape[0])
        assert_equal(cg['stars', 'particle_mass'].shape[0],
                     grid['stars', 'particle_mass'].shape[0])
