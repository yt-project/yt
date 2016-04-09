import yt
from yt.analysis_modules.cosmological_observation.api import \
     LightCone

# Create a LightCone object extending from z = 0 to z = 0.1.

# We have already set up the redshift dumps to be
# used for this, so we will not use any of the time
# data dumps.
lc = LightCone('enzo_tiny_cosmology/32Mpc_32.enzo',
               'Enzo', 0., 0.1,
               observer_redshift=0.0,
               time_data=False)

# Calculate a randomization of the solution.
lc.calculate_light_cone_solution(seed=123456789, filename="LC/solution.txt")

# Choose the field to be projected.
field = 'szy'

# Use the LightCone object to make a projection with a 600 arcminute
# field of view and a resolution of 60 arcseconds.
# Set njobs to -1 to have one core work on each projection
# in parallel.
lc.project_light_cone((600.0, "arcmin"), (60.0, "arcsec"), field,
                      weight_field=None,
                      save_stack=True,
                      save_final_image=True,
                      save_slice_images=True,
                      njobs=-1)

# By default, the light cone projections are kept in the LC directory,
# but this moves them back to the current directory so that they're rendered
# in our cookbook.
import shutil, glob
for file in glob.glob('LC/*png'):
    shutil.move(file, '.')
