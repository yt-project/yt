from yt.mods import *
from yt.analysis_modules.cosmological_observation.light_cone.api import *

# Instantiate a light cone.
lc = LightCone("enzo_tiny_cosmology/32Mpc_32.enzo", 'Enzo', 0, 0.1,
               observer_redshift=0.0,
               field_of_view_in_arcminutes=120.0,
               image_resolution_in_arcseconds=60.0,
               use_minimum_datasets=True,
               time_data=False,
               output_dir='LC_U', output_prefix='LightCone')

# Try to find 10 solutions that have at most 10% volume in
# common and give up after 50 consecutive failed attempts.
# The recycle=True setting tells the code to first attempt
# to use solutions with the same projection axes as other
# solutions.  This will save time when making the projection.
find_unique_solutions(lc, max_overlap=0.10, failures=50,
                      seed=123456789, recycle=True,
                      solutions=10, filename='LC_U/unique.dat')

# Choose the field to be projected.
field = 'SZY'

# Make light cone projections with each of the random seeds
# found above.  All output files will be written with unique
# names based on the random seed numbers.
project_unique_light_cones(lc, 'LC_U/unique.dat', field,
                           save_stack=False,
                           save_final_image=True,
                           save_slice_images=False)
