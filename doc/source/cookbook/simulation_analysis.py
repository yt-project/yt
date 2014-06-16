from yt.mods import *

# Instantiate a time series object for an Enzo simulation..
my_sim = simulation('enzo_tiny_cosmology/32Mpc_32.enzo', 'Enzo')

# Get a time series for all data made by the simulation.
my_sim.get_time_series()

# Calculate and store extrema for all datasets.
all_storage = {}
for my_storage, ds in my_sim.piter(storage=all_storage):
    all_data = ds.all_data()
    my_extrema = all_data.quantities['Extrema']('density')

    # Save to storage so we can get at it afterward.
    my_storage.result = my_extrema

# Print out all the values we calculated.
for my_result in all_storage.values():
    print my_result
