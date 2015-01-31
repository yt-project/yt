import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog

# Load the data set with the full simulation information
# and rockstar halos
data_ds = yt.load('Enzo_64/RD0006/RedshiftOutput0006')
halos_ds = yt.load('rockstar_halos/halos_0.0.bin')

# Instantiate a catalog using those two paramter files
hc = HaloCatalog(data_ds=data_ds, halos_ds=halos_ds)

# Filter out less massive halos
hc.add_filter("quantity_value", "particle_mass", ">", 1e14, "Msun")

# attach a sphere object to each halo whose radius extends
#   to twice the radius of the halo
hc.add_callback("sphere", factor=2.0)

# use the sphere to calculate radial profiles of gas density
# weighted by cell volume in terms of the virial radius
hc.add_callback("profile", ["radius"],
                [("gas", "overdensity")],
                weight_field="cell_volume",
                accumulation=True,
                storage="virial_quantities_profiles")


hc.add_callback("virial_quantities", ["radius"],
                profile_storage="virial_quantities_profiles")
hc.add_callback('delete_attribute', 'virial_quantities_profiles')

field_params = dict(virial_radius=('quantity', 'radius_200'))
hc.add_callback('sphere', radius_field='radius_200', factor=5,
                field_parameters=field_params)
hc.add_callback('profile', ['virial_radius_fraction'], 
                [('gas', 'temperature')],
                storage='virial_profiles',
                weight_field='cell_mass',
                accumulation=False, output_dir='profiles')

# Save the profiles
hc.add_callback("save_profiles", storage="virial_profiles",
                output_dir="profiles")

hc.create()
