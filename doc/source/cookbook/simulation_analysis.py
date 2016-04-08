import yt
yt.enable_parallelism()
import collections

# Enable parallelism in the script (assuming it was called with
# `mpirun -np <n_procs>` )
yt.enable_parallelism()

# By using wildcards such as ? and * with the load command, we can load up a
# Time Series containing all of these datasets simultaneously.
ts = yt.load('enzo_tiny_cosmology/DD????/DD????')

# Calculate and store density extrema for all datasets along with redshift
# in a data dictionary with entries as tuples

# Create an empty dictionary
data = {}

# Iterate through each dataset in the Time Series (using piter allows it
# to happen in parallel automatically across available processors)
for ds in ts.piter():
    ad = ds.all_data()
    extrema = ad.quantities.extrema('density')

    # Fill the dictionary with extrema and redshift information for each dataset
    data[ds.basename] = (extrema, ds.current_redshift)

# Convert dictionary to ordered dictionary to get the right order
od = collections.OrderedDict(sorted(data.items()))

# Print out all the values we calculated.
print("Dataset      Redshift        Density Min      Density Max")
print("---------------------------------------------------------")
for key, val in od.items():
    print("%s       %05.3f          %5.3g g/cm^3   %5.3g g/cm^3" % \
           (key, val[1], val[0][0], val[0][1]))
