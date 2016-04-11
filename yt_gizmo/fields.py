"""
Gizmo-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Britton Smith.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.species_fields import \
    add_species_field_by_fraction
from yt.frontends.gadget.fields import \
    GadgetFieldInfo

class GizmoFieldInfo(GadgetFieldInfo):
    def setup_gas_particle_fields(self, ptype):
        super(GizmoFieldInfo, self).setup_gas_particle_fields(ptype)
        self.alias((ptype, "temperature"), (ptype, "Temperature"))

        self.alias((ptype, "H_fraction"), (ptype, "NeutralHydrogenAbundance"))
        add_species_field_by_fraction(self, ptype, "H", particle_type=True)
