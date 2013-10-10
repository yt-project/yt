"""
OWLS-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import *
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    NullFunc, \
    TranslationFunc
import yt.fields.universal_fields
from .definitions import \
    gadget_ptypes, \
    ghdf5_ptypes
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_scalar_functions, \
    _field_concat, _field_concat_slice
from yt.utilities.physical_constants import \
    mass_sun_cgs

OWLSFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_owls_field = OWLSFieldInfo.add_field

KnownOWLSFields = FieldInfoContainer()
add_owls_field = KnownOWLSFields.add_field

GadgetFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_gadget_field = GadgetFieldInfo.add_field

KnownGadgetFields = FieldInfoContainer()
add_gadget_field = KnownGadgetFields.add_field

GadgetHDF5FieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
KnownGadgetHDF5Fields = FieldInfoContainer()

TipsyFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_Tipsy_field = TipsyFieldInfo.add_field

KnownTipsyFields = FieldInfoContainer()
add_tipsy_field = KnownTipsyFields.add_field

# Here are helper functions for things like vector fields and so on.

def _get_conv(cf):
    def _convert(data):
        return data.convert(cf)
    return _convert

# TIPSY
# =====

def _setup_particle_fields(registry, ptype):
    registry.add_field((ptype, "Mass"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("mass"),
        units = r"\mathrm{g}")
    registry.add_field((ptype, "Velocities"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("velocity"),
        units = r"\mathrm{cm}/\mathrm{s}")
    # Note that we have to do this last so that TranslationFunc operates
    # correctly.
    particle_deposition_functions(ptype, "Coordinates", "Mass",
                                  registry)
    particle_scalar_functions(ptype, "Coordinates", "Velocities",
                              registry)
    for fname in ["Coordinates", "ParticleIDs", "Epsilon", "Phi"]:
        registry.add_field((ptype, fname), function=NullFunc,
                particle_type = True)

def SmoothedGas(field, data):
    pos = data["PartType0", "Coordinates"]
    sml = data["PartType0", "SmoothingLength"]
    dens = data["PartType0", "Density"]
    rv = data.deposit(pos, [sml, dens], method="simple_smooth")
    return rv
OWLSFieldInfo.add_field(("deposit", "PartType0_simple_smooth"),
                function = SmoothedGas, validators = [ValidateSpatial()])

