import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")  # load data

field = "temperature"  # The field to average
weight = "cell_mass"  # The weight for the average

dd = ds.all_data()  # This is a region describing the entire box,
                      # but note it doesn't read anything in yet!
# We now use our 'quantities' call to get the average quantity
average_value = dd.quantities["WeightedAverageQuantity"](field, weight)

print "Average %s (weighted by %s) is %0.5e" % (field, weight, average_value)
