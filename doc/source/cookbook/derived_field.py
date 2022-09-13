import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# You can create a derived field by manipulating any existing derived fields
# in any way you choose.  In this case, let's just make a simple one:
# thermal_energy_density = 3/2 nkT

# First create a function which yields your new derived field
def thermal_energy_dens(field, data):
    return (3 / 2) * data["gas", "number_density"] * data["gas", "kT"]


# Then add it to your dataset and define the units
ds.add_field(
    ("gas", "thermal_energy_density"),
    units="erg/cm**3",
    function=thermal_energy_dens,
    sampling_type="cell",
)

# It will now show up in your derived_field_list
for i in sorted(ds.derived_field_list):
    print(i)

# Let's use it to make a projection
ad = ds.all_data()
yt.ProjectionPlot(
    ds,
    "x",
    ("gas", "thermal_energy_density"),
    weight_field=("gas", "density"),
    width=(200, "kpc"),
).save()
