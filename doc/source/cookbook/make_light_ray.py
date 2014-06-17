import os
import sys
import yt
from yt.analysis_modules.halo_profiler.api import HaloProfiler
from yt.analysis_modules.cosmological_observation.light_ray.light_ray import \
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

# Configure the HaloProfiler.
# These are keyword arguments given when creating a
# HaloProfiler object.
halo_profiler_kwargs = {'halo_list_file': 'HopAnalysis.out',
                        'output_dir' : 'halo_analysis'}

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

# Make a light ray, and set njobs to -1 to use one core
# per dataset.
lr.make_light_ray(seed=123456789,
                  solution_filename='LR/lightraysolution.txt',
                  data_filename='LR/lightray.h5',
                  fields=['temperature', 'density'],
                  get_nearest_halo=True,
                  nearest_halo_fields=['TotalMassMsun_100',
                                       'RadiusMpc_100'],
                  halo_profiler_parameters=halo_profiler_parameters,
                  get_los_velocity=True,
                  njobs=-1)
