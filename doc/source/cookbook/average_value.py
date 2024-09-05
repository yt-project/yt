import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")  # load data

field = ("gas", "temperature")  # The field to average
weight = ("gas", "mass")  # The weight for the average

ad = ds.all_data()  # This is a region describing the entire box,
# but note it doesn't read anything in yet!

# We now use our 'quantities' call to get the average quantity
average_value = ad.quantities.weighted_average_quantity(field, weight)

print(
    "Average %s (weighted by %s) is %0.3e %s"
    % (field[1], weight[1], average_value, average_value.units)
)
