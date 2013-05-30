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
    ValidateGridType
import yt.data_objects.universal_fields

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
             projection_conversion = '1')

    def particle_density(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        d /= data["CellVolume"]
        return d

    registry.add_field(("deposit", "%s_density" % ptype),
             function = particle_density,
             validators = [ValidateSpatial()],
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projection_conversion = 'cm')

for ptype in ["Gas", "DarkMatter", "Stars"]:
    _particle_functions(ptype, "Coordinates", "Mass", TipsyFieldInfo)
   
# GADGET
# ======

# Among other things we need to set up Coordinates

_gadget_ptypes = ("Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry")

def _gadget_particle_fields(ptype):
    def _Mass(field, data):
        pind = _gadget_ptypes.index(ptype)
        if data.pf["Massarr"][pind] == 0.0:
            return data[ptype, "Masses"]
        mass = np.ones(data[ptype, "Coordinates"].shape[0], dtype="float64")
        mass *= data.pf["Massarr"][pind]
        return mass
    GadgetFieldInfo.add_field((ptype, "Mass"), function=_Mass,
                              particle_type = True)

for ptype in _gadget_ptypes:
    _gadget_particle_fields(ptype)
    _particle_functions(ptype, "Coordinates", "Mass", GadgetFieldInfo)
