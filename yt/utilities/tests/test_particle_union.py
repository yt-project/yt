from __future__ import print_function
import yt
from yt.testing import requires_file
from yt.data_objects.particle_filters import add_particle_filter


# Dataset required for this test
iso_galaxy = 'IsolatedGalaxy/galaxy0030/galaxy0030'


# Stars function for stars particle filter
# Define particle filter to filter stars younger than 100 Myr
def Stars(pfilter, data):
    return (data.ds.current_time - data["all",
                                        "creation_time"]).in_units('Myr') > 0

# Simple test to create a particle filter, then access the deposition field
# (was a problem in previous versions on this dataset because there are chunks
#  with no stars in them).


@requires_file(iso_galaxy)
def test_particle_filter():
    add_particle_filter("stars", function=Stars, filtered_type='all',
                        requires=["creation_time"])
    ds = yt.load(iso_galaxy)
    ds.add_particle_filter('stars')
    ad = ds.all_data()
    print(ad['deposit', 'stars_cic'])
