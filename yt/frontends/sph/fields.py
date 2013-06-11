"""
OWLS-specific fields

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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
import yt.data_objects.universal_fields
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

TipsyFieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_Tipsy_field = TipsyFieldInfo.add_field

KnownTipsyFields = FieldInfoContainer()
add_tipsy_field = KnownTipsyFields.add_field

def _particle_functions(ptype, coord_name, mass_name, registry):
    def particle_count(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, method = "count")
        return d
    registry.add_field(("deposit", "%s_count" % ptype),
             function = particle_count,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Count}" % ptype,
             projection_conversion = '1')

    def particle_mass(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        return d

    registry.add_field(("deposit", "%s_mass" % ptype),
             function = particle_mass,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Mass}" % ptype,
             units = r"\mathrm{g}",
             projected_units = r"\mathrm{g}\/\mathrm{cm}",
             projection_conversion = 'cm')

    def particle_density(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        d /= data["CellVolume"]
        return d

    registry.add_field(("deposit", "%s_density" % ptype),
             function = particle_density,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    def particle_cic(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
        d /= data["CellVolume"]
        return d

    registry.add_field(("deposit", "%s_cic" % ptype),
             function = particle_cic,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s CIC Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    # Now some translation functions.

    registry.add_field((ptype, "ParticleMass"),
            function = TranslationFunc((ptype, mass_name)),
            particle_type = True,
            units = r"\mathrm{g}")

    def _ParticleMassMsun(field, data):
        return data[ptype, mass_name].copy()
    def _conv_Msun(data):
        return 1.0/mass_sun_cgs

    registry.add_field((ptype, "ParticleMassMsun"),
            function = _ParticleMassMsun,
            convert_function = _conv_Msun,
            particle_type = True,
            units = r"\mathrm{M}_\odot")

    # For 'all', which is a special field, we skip adding a few types.
    
    if ptype == "all": return

    # Now we have to set up the various velocity and coordinate things.  In the
    # future, we'll actually invert this and use the 3-component items
    # elsewhere, and stop using these.
    
    # Note that we pass in _ptype here so that it's defined inside the closure.
    def _get_coord_funcs(axi, _ptype):
        def _particle_velocity(field, data):
            return data[_ptype, "Velocities"][:,axi]
        def _particle_position(field, data):
            return data[_ptype, "Coordinates"][:,axi]
        return _particle_velocity, _particle_position
    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs(axi, ptype)
        registry.add_field((ptype, "particle_velocity_%s" % ax),
            particle_type = True, function = v)
        registry.add_field((ptype, "particle_position_%s" % ax),
            particle_type = True, function = p)

# Here are helper functions for things like vector fields and so on.

def _get_conv(cf):
    def _convert(data):
        return data.convert(cf)
    return _convert

def _field_concat(fname):
    def _AllFields(field, data):
        v = []
        for ptype in data.pf.particle_types:
            if ptype == "all": continue
            v.append(data[ptype, fname].copy())
        rv = np.concatenate(v, axis=0)
        return rv
    return _AllFields

def _field_concat_slice(fname, axi):
    def _AllFields(field, data):
        v = []
        for ptype in data.pf.particle_types:
            if ptype == "all": continue
            v.append(data[ptype, fname][:,axi])
        rv = np.concatenate(v, axis=0)
        return rv
    return _AllFields

# TIPSY
# =====

for ptype in ["Gas", "DarkMatter", "Stars"]:
    KnownTipsyFields.add_field((ptype, "Mass"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("mass"),
        units = r"\mathrm{g}")
    KnownTipsyFields.add_field((ptype, "Velocities"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("velocity"),
        units = r"\mathrm{cm}/\mathrm{s}")
    # Note that we have to do this last so that TranslationFunc operates
    # correctly.
    _particle_functions(ptype, "Coordinates", "Mass", TipsyFieldInfo)
_particle_functions("all", "Coordinates", "Mass", TipsyFieldInfo)

for fname in ["Coordinates", "Velocities", "ParticleIDs", "Mass",
              "Epsilon", "Phi"]:
    func = _field_concat(fname)
    TipsyFieldInfo.add_field(("all", fname), function=func,
            particle_type = True)

for iname, oname in [("Coordinates", "particle_position_"),
                     ("Velocities", "particle_velocity_")]:
    for axi, ax in enumerate("xyz"):
        func = _field_concat_slice(iname, axi)
        TipsyFieldInfo.add_field(("all", oname + ax), function=func,
                particle_type = True)

# GADGET
# ======

# Among other things we need to set up Coordinates

_gadget_ptypes = ("Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry")

# This has to be done manually for Gadget, because some of the particles will
# have uniform mass
def _gadget_particle_fields(ptype):
    def _Mass(field, data):
        pind = _gadget_ptypes.index(ptype)
        if data.pf["Massarr"][pind] == 0.0:
            return data[ptype, "Masses"].copy()
        mass = np.ones(data[ptype, "ParticleIDs"].shape[0], dtype="float64")
        # Note that this is an alias, which is why we need to apply conversion
        # here.  Otherwise we'd have an asymmetry.
        mass *= data.pf["Massarr"][pind] * data.convert("mass")
        return mass
    GadgetFieldInfo.add_field((ptype, "Mass"), function=_Mass,
                              particle_type = True)

for fname in ["Coordinates", "Velocities", "ParticleIDs",
              # Note: Mass, not Masses
              "Mass"]:
    func = _field_concat(fname)
    GadgetFieldInfo.add_field(("all", fname), function=func,
            particle_type = True)

for ptype in _gadget_ptypes:
    KnownGadgetFields.add_field((ptype, "Masses"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("mass"),
        units = r"\mathrm{g}")
    _gadget_particle_fields(ptype)
    KnownGadgetFields.add_field((ptype, "Velocities"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("velocity"),
        units = r"\mathrm{cm}/\mathrm{s}")
    _particle_functions(ptype, "Coordinates", "Mass", GadgetFieldInfo)
    KnownGadgetFields.add_field((ptype, "Coordinates"), function=NullFunc,
        particle_type = True)
_particle_functions("all", "Coordinates", "Mass", GadgetFieldInfo)

# Now we have to manually apply the splits for "all", since we don't want to
# use the splits defined above.

for iname, oname in [("Coordinates", "particle_position_"),
                     ("Velocities", "particle_velocity_")]:
    for axi, ax in enumerate("xyz"):
        func = _field_concat_slice(iname, axi)
        GadgetFieldInfo.add_field(("all", oname + ax), function=func,
                particle_type = True)

# OWLS
# ====

# I am optimistic that some day we will be able to get rid of much of this, and
# make OWLS a subclass of Gadget fields.

_owls_ptypes = ("PartType0", "PartType1", "PartType2", "PartType3",
                "PartType4")

for fname in ["Coordinates", "Velocities", "ParticleIDs",
              # Note: Mass, not Masses
              "Mass"]:
    func = _field_concat(fname)
    OWLSFieldInfo.add_field(("all", fname), function=func,
            particle_type = True)

def _owls_particle_fields(ptype):
    def _Mass(field, data):
        pind = _owls_ptypes.index(ptype)
        if data.pf["MassTable"][pind] == 0.0:
            raise RuntimeError
        mass = np.ones(data[ptype, "ParticleIDs"].shape[0], dtype="float64")
        # Note that this is an alias, which is why we need to apply conversion
        # here.  Otherwise we'd have an asymmetry.
        mass *= data.pf["MassTable"][pind] 
        return mass
    OWLSFieldInfo.add_field((ptype, "Mass"), function=_Mass,
                            convert_function = _get_conv("mass"),
                            particle_type = True)

for ptype in _owls_ptypes:
    # Note that this adds a "Known" Mass field and a "Derived" Mass field.
    # This way the "Known" will get used, and if it's not there, it will use
    # the derived.
    KnownOWLSFields.add_field((ptype, "Mass"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("mass"),
        units = r"\mathrm{g}")
    _owls_particle_fields(ptype)
    KnownOWLSFields.add_field((ptype, "Velocities"), function=NullFunc,
        particle_type = True,
        convert_function=_get_conv("velocity"),
        units = r"\mathrm{cm}/\mathrm{s}")
    _particle_functions(ptype, "Coordinates", "Mass", OWLSFieldInfo)
    KnownOWLSFields.add_field((ptype, "Coordinates"), function=NullFunc,
        particle_type = True)
_particle_functions("all", "Coordinates", "Mass", OWLSFieldInfo)

# Now we have to manually apply the splits for "all", since we don't want to
# use the splits defined above.

for iname, oname in [("Coordinates", "particle_position_"),
                     ("Velocities", "particle_velocity_")]:
    for axi, ax in enumerate("xyz"):
        func = _field_concat_slice(iname, axi)
        OWLSFieldInfo.add_field(("all", oname + ax), function=func,
                particle_type = True)
