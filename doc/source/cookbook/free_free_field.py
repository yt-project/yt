### THIS RECIPE IS CURRENTLY BROKEN IN YT-3.0
### DO NOT TRUST THIS RECIPE UNTIL THIS LINE IS REMOVED

import numpy as np
import yt
# Need to grab the proton mass from the constants database
from yt.utilities.physical_constants import mp

exit()
# Define the emission field

keVtoerg = 1.602e-9  # Convert energy in keV to energy in erg
KtokeV = 8.617e-08  # Convert degrees Kelvin to degrees keV
sqrt3 = np.sqrt(3.)
expgamma = 1.78107241799  # Exponential of Euler's constant


def _FreeFree_Emission(field, data):

    if data.has_field_parameter("Z"):
        Z = data.get_field_parameter("Z")
    else:
        Z = 1.077  # Primordial H/He plasma

    if data.has_field_parameter("mue"):
        mue = data.get_field_parameter("mue")
    else:
        mue = 1./0.875  # Primordial H/He plasma

    if data.has_field_parameter("mui"):
        mui = data.get_field_parameter("mui")
    else:
        mui = 1./0.8125  # Primordial H/He plasma

    if data.has_field_parameter("Ephoton"):
        Ephoton = data.get_field_parameter("Ephoton")
    else:
        Ephoton = 1.0  # in keV

    if data.has_field_parameter("photon_emission"):
        photon_emission = data.get_field_parameter("photon_emission")
    else:
        photon_emission = False  # Flag for energy or photon emission

    n_e = data["density"]/(mue*mp)
    n_i = data["density"]/(mui*mp)
    kT = data["temperature"]*KtokeV

    # Compute the Gaunt factor

    g_ff = np.zeros(kT.shape)
    g_ff[Ephoton/kT > 1.] = np.sqrt((3./np.pi)*kT[Ephoton/kT > 1.]/Ephoton)
    g_ff[Ephoton/kT < 1.] = (sqrt3/np.pi)*np.log((4./expgamma) *
                                                 kT[Ephoton/kT < 1.]/Ephoton)

    eps_E = 1.64e-20*Z*Z*n_e*n_i/np.sqrt(data["temperature"]) * \
        np.exp(-Ephoton/kT)*g_ff

    if photon_emission:
        eps_E /= (Ephoton*keVtoerg)

    return eps_E

yt.add_field("FreeFree_Emission", function=_FreeFree_Emission)

# Define the luminosity derived quantity
def _FreeFreeLuminosity(data):
    return (data["FreeFree_Emission"]*data["cell_volume"]).sum()


def _combFreeFreeLuminosity(data, luminosity):
    return luminosity.sum()

yt.add_quantity("FreeFree_Luminosity", function=_FreeFreeLuminosity,
                combine_function=_combFreeFreeLuminosity, n_ret=1)

ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150")

sphere = ds.sphere(ds.domain_center, (100., "kpc"))

# Print out the total luminosity at 1 keV for the sphere

print("L_E (1 keV, primordial) = ", sphere.quantities["FreeFree_Luminosity"]())

# The defaults for the field assume a H/He primordial plasma.
# Let's set the appropriate parameters for a pure hydrogen plasma.

sphere.set_field_parameter("mue", 1.0)
sphere.set_field_parameter("mui", 1.0)
sphere.set_field_parameter("Z", 1.0)

print("L_E (1 keV, pure hydrogen) = ", sphere.quantities["FreeFree_Luminosity"]())

# Now let's print the luminosity at an energy of E = 10 keV

sphere.set_field_parameter("Ephoton", 10.0)

print("L_E (10 keV, pure hydrogen) = ", sphere.quantities["FreeFree_Luminosity"]())

# Finally, let's set the flag for photon emission, to get the total number
# of photons emitted at this energy:

sphere.set_field_parameter("photon_emission", True)

print("L_ph (10 keV, pure hydrogen) = ", sphere.quantities["FreeFree_Luminosity"]())
