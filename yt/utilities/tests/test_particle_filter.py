import numpy as np
import yt
from yt.mods import *
from yt.testing import *
from yt.utilities.particle_generator import *
from yt.frontends.stream.api import load_uniform_grid, refine_amr
import yt.utilities.initial_conditions as ic
import yt.utilities.flagging_methods as fm
from IPython import embed
from yt.units.yt_array import uconcatenate
from yt.data_objects.particle_filters import add_particle_filter

def setup() :
    pass

# Stars function for stars particle filter
# Define particle filter to filter stars younger than 100 Myr
def Stars(pfilter, data):
    filter = (data.ds.current_time - data["all", "creation_time"]).in_units('Myr') > 0
    return filter

# Simple test to create a particle filter, then access the deposition field 
# (was a problem in previous versions on this dataset because there are chunks 
#  with no stars in them).
def test_particle_filter() :
    add_particle_filter("stars", function=Stars, filtered_type='all', requires=["creation_time"])
    ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
    ds.add_particle_filter('stars')
    ad = ds.all_data()
    print ad['deposit', 'stars_cic']
