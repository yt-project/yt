import os
import yt
from yt.analysis_modules.cosmological_observation.api import \
    LightRay

# Create a directory for the light rays
if not os.path.isdir("LR"): 
    os.mkdir('LR')
     
# Create a LightRay object extending from z = 0 to z = 0.1
# and use only the redshift dumps.
lr = LightRay("enzo_tiny_cosmology/32Mpc_32.enzo",
              'Enzo', 0.0, 0.1,
              use_minimum_datasets=True,
              time_data=False)

# Make a light ray, and set njobs to -1 to use one core
# per dataset.
lr.make_light_ray(seed=123456789,
                  solution_filename='LR/lightraysolution.txt',
                  data_filename='LR/lightray.h5',
                  fields=['temperature', 'density'],
                  get_los_velocity=True,
                  njobs=-1)
