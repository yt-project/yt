import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# You can create a derived field by manipulating any existing derived fields
# in any way you choose.  In this case, let's just make a simple one:
# number_density = density / mass

# First create a function which yields your new derived field
def number_density(field, data):
    return data['gas', 'density']/data['gas', 'mass']

# Then add it to your dataset and define the units
ds.add_field(("gas", "number_density"), units="cm**-3", function=number_density)

# It will now show up in your derived_field_list
for i in sorted(ds.derived_field_list): 
    print i

# Let's use it to make a projection of the entire volume!
ad = ds.all_data()
yt.ProjectionPlot(ds, "x", "number_density").save()
