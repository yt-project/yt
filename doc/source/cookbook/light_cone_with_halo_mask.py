### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import yt

from yt.analysis_modules.cosmological_observation.light_cone.light_cone import LightCone
from yt.analysis_modules.halo_profiler.api import HaloProfiler

# Instantiate a light cone object as usual.
lc = LightCone('enzo_tiny_cosmology/32Mpc_32.enzo',
               'Enzo', 0, 0.1,
               observer_redshift=0.0,
               field_of_view_in_arcminutes=600.0,
               image_resolution_in_arcseconds=60.0,
               time_data=False,
               output_dir='LC_HM', output_prefix='LightCone')

# Calculate the light cone solution.
lc.calculate_light_cone_solution(seed=123456789,
                                 filename='LC_HM/lightcone.dat')


# Configure the HaloProfiler.
# These are keyword arguments given when creating a
# HaloProfiler object.
halo_profiler_kwargs = {'halo_list_file': 'HopAnalysis.out',
                        'output_dir': 'halo_analysis'}

# Create a list of actions for the HaloProfiler to take.
halo_profiler_actions = []

# Each item in the list is a dictionary containing three things:
# 1. 'function' - the function to be called.
# 2. 'args' - a list of arguments given with the function.
# 3. 'kwargs' - a dictionary of keyword arguments.

# Add a virial filter.
halo_profiler_actions.append({'function': HaloProfiler.add_halo_filter,
                              'args': [VirialFilter],
                              'kwargs': {'must_be_virialized':False,
                                         'overdensity_field':'ActualOverdensity',
                                         'virial_overdensity':100,
                                         'virial_filters':[['TotalMassMsun','>','1e5']],
                                         'virial_quantities':['TotalMassMsun','RadiusMpc']}})

# Add a call to make the profiles.
halo_profiler_actions.append({'function': HaloProfiler.make_profiles,
                              'kwargs': {'filename': "VirializedHalos.out"}})

# Specify the desired halo list is the filtered list.
# If 'all' is given instead, the full list will be used.
halo_list = 'filtered'

# Put them all into one dictionary.
halo_profiler_parameters=dict(halo_profiler_kwargs=halo_profiler_kwargs,
                              halo_profiler_actions=halo_profiler_actions,
                              halo_list=halo_list)

# Get the halo list for the active solution of this light cone using
# the HaloProfiler settings set up above.
# Write the boolean map to an hdf5 file called 'halo_mask.h5'.
# Write a text file detailing the location, redshift, radius, and mass
# of each halo in light cone projection.
lc.get_halo_mask(mask_file='LC_HM/halo_mask.h5',
                 map_file='LC_HM/halo_map.out',
                 cube_file='LC_HM/halo_cube.h5',
                 virial_overdensity=100,
                 halo_profiler_parameters=halo_profiler_parameters,
                 njobs=1, dynamic=False)

# Choose the field to be projected.
field = 'SZY'

# Make the light cone projection and apply the halo mask.
pc = lc.project_light_cone(field, save_stack=False,
                           save_final_image=True,
                           save_slice_images=False,
                           apply_halo_mask=True)
