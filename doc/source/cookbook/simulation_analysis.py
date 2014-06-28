import yt
import collections

# Instantiate a time series object for an Enzo simulation..
sim = yt.simulation('enzo_tiny_cosmology/32Mpc_32.enzo', 'Enzo')

# Get a time series for all data made by the simulation.
sim.get_time_series()

# Calculate and store extrema for all datasets along with redshift
# in a data dictionary with entries as tuples

# Note that by using sim.piter(), we are automatically 
# forcing yt to do this in parallel
data = {}
for ds in sim.piter():
    ad = ds.all_data()
    extrema = ad.quantities.extrema('density')
    data[ds.basename] = (extrema, ds.current_redshift)

# Convert dictionary to ordered dictionary to get the right order
od = collections.OrderedDict(sorted(data.items()))

# Print out all the values we calculated.
print "Dataset      Redshift        Density Min      Density Max"
print "---------------------------------------------------------"
for k, v in od.iteritems(): 
    print "%s       %05.3f          %5.3g g/cm^3   %5.3g g/cm^3" % (k, v[1], v[0][0], v[0][1])
