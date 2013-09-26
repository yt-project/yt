"""
RAMSES-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os

from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    NullFunc, \
    TranslationFunc, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
from yt.data_objects.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs, \
    mass_sun_cgs, \
    mh
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator
import yt.utilities.fortran_utils as fpu
from yt.funcs import mylog
import numpy as np

RAMSESFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo, "RFI")
add_field = RAMSESFieldInfo.add_field

KnownRAMSESFields = FieldInfoContainer()
add_ramses_field = KnownRAMSESFields.add_field

known_ramses_fields = [
    "Density",
    "x-velocity",
    "y-velocity",
    "z-velocity",
    "Pressure",
    "Metallicity",
]

for f in known_ramses_fields:
    if f not in KnownRAMSESFields:
        add_ramses_field(f, function=NullFunc, take_log=True,
                  validators = [ValidateDataField(f)])

def dx(field, data):
    return data.fwidth[:,0]
add_field("dx", function=dx)

def dy(field, data):
    return data.fwidth[:,1]
add_field("dy", function=dy)

def dz(field, data):
    return data.fwidth[:,2]
add_field("dz", function=dz)

def _convertDensity(data):
    return data.convert("Density")
KnownRAMSESFields["Density"]._units = r"\rm{g}/\rm{cm}^3"
KnownRAMSESFields["Density"]._projected_units = r"\rm{g}/\rm{cm}^2"
KnownRAMSESFields["Density"]._convert_function=_convertDensity

def _convertPressure(data):
    return data.convert("Pressure")
KnownRAMSESFields["Pressure"]._units=r"\rm{dyne}/\rm{cm}^{2}/\mu"
KnownRAMSESFields["Pressure"]._convert_function=_convertPressure

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = KnownRAMSESFields["%s-velocity" % ax]
    f._units = r"\rm{cm}/\rm{s}"
    f._convert_function = _convertVelocity
    f.take_log = False

known_ramses_particle_fields = [
    "particle_position_x",
    "particle_position_y",
    "particle_position_z",
    "particle_velocity_x",
    "particle_velocity_y",
    "particle_velocity_z",
    "particle_mass",
    "particle_identifier",
    "particle_refinement_level",
    "particle_age",
    "particle_metallicity",
]

for f in known_ramses_particle_fields:
    add_ramses_field(("all", f), function=NullFunc, take_log=True,
              particle_type = True)

for ax in 'xyz':
    KnownRAMSESFields["all", "particle_velocity_%s" % ax]._convert_function = \
        _convertVelocity

def _convertParticleMass(data):
    return data.convert("mass")

KnownRAMSESFields["all", "particle_mass"]._convert_function = \
        _convertParticleMass
KnownRAMSESFields["all", "particle_mass"]._units = r"\mathrm{g}"

def _Temperature(field, data):
    rv = data["Pressure"]/data["Density"]
    rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
    return rv
add_field("Temperature", function=_Temperature, units=r"\rm{K}")

# We'll add a bunch of species fields here.  In the not too distant future,
# we'll be moving all of these to a unified field location, so they can be
# shared between various frontends.

# NOTE: No Electron here because I don't know how RAMSES handles them, and if
# they are handled differently than Enzo does (where they are scaled to mh)

_speciesList = ["HI", "HII",
                "HeI", "HeII", "HeIII",
                "H2I", "H2II", "HM",
                "DI", "DII", "HDI"]
_speciesMass = {"HI": 1.0, "HII": 1.0,
                "HeI": 4.0, "HeII": 4.0, "HeIII": 4.0,
                "H2I": 2.0, "H2II": 2.0, "HM": 1.0,
                "DI": 2.0, "DII": 2.0, "HDI": 3.0}

def _SpeciesComovingDensity(field, data):
    sp = field.name.split("_")[0] + "_Density"
    ef = (1.0 + data.pf.current_redshift)**3.0
    return data[sp] / ef

def _SpeciesFraction(field, data):
    species = field.name.split("_")[0]
    sp = "%s_NumberDensity" % species
    return mh * data[sp] * _speciesMass[species]

def _SpeciesMass(field, data):
    sp = field.name.split("_")[0] + "_Density"
    return data[sp] * data["CellVolume"]

def _SpeciesDensity(field, data):
    species = field.name.split("_")[0]
    sp = "%s_NumberDensity" % species
    return mh * data[sp] * _speciesMass[species] * data["Density"]

def _convertCellMassMsun(data):
    return 1.0/mass_sun_cgs # g^-1

for species in _speciesList:
    add_field("%s_Fraction" % species,
             function = _SpeciesFraction,
             display_name = "%s\/Fraction" % species)
    add_field("%s_Density" % species,
             function = _SpeciesDensity,
             display_name = "%s\/Density" % species,
             units = r"\rm{g}/\rm{cm}^3",
             projected_units = r"\rm{g}/\rm{cm}^2")
    add_field("Comoving_%s_Density" % species,
             function=_SpeciesComovingDensity,
             validators=ValidateDataField("%s_Density" % species),
             display_name="Comoving\/%s\/Density" % species)
    add_field("%s_Mass" % species, units=r"\rm{g}", 
              function=_SpeciesMass, 
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)
    add_field("%s_MassMsun" % species, units=r"M_{\odot}", 
              function=_SpeciesMass, 
              convert_function=_convertCellMassMsun,
              validators=ValidateDataField("%s_Density" % species),
              display_name="%s\/Mass" % species)

# PARTICLE FIELDS
particle_vector_functions("all", ["particle_position_%s" % ax for ax in 'xyz'],
                                 ["particle_velocity_%s" % ax for ax in 'xyz'],
                          RAMSESFieldInfo)
particle_deposition_functions("all", "Coordinates", "particle_mass",
                               RAMSESFieldInfo)
_cool_axes = ("lognH", "logT", "logTeq")
_cool_arrs = ("metal", "cool", "heat", "metal_prime", "cool_prime",
              "heat_prime", "mu", "abundances")
_cool_species = ("Electron_NumberDensity",
                 "HI_NumberDensity",
                 "HII_NumberDensity",
                 "HeI_NumberDensity",
                 "HeII_NumberDensity",
                 "HeIII_NumberDensity")

def create_cooling_fields(filename, field_info):
    if not os.path.exists(filename): return
    def _create_field(name, interp_object):
        def _func(field, data):
            shape = data["Temperature"].shape
            d = {'lognH': np.log10(data["Density"]/mh).ravel(),
                 'logT' : np.log10(data["Temperature"]).ravel()}
            rv = 10**interp_object(d).reshape(shape)
            return rv
        field_info.add_field(name = name, function=_func,
                             units = r"\rm{cm}^{-3}",
                             projected_units = r"\rm{cm}^{-2}")
    avals = {}
    tvals = {}
    with open(filename, "rb") as f:
        n1, n2 = fpu.read_vector(f, 'i')
        n = n1 * n2
        for ax in _cool_axes:
            avals[ax] = fpu.read_vector(f, 'd')
        for tname in _cool_arrs:
            var = fpu.read_vector(f, 'd')
            if var.size == n1*n2:
                tvals[tname] = var.reshape((n1, n2), order='F')
            else:
                var = var.reshape((n1, n2, var.size / (n1*n2)), order='F')
                for i in range(var.shape[-1]):
                    tvals[_cool_species[i]] = var[:,:,i]
    
    for n in tvals:
        interp = BilinearFieldInterpolator(tvals[n],
                    (avals["lognH"], avals["logT"]),
                    ["lognH", "logT"], truncate = True)
        _create_field(n, interp)
