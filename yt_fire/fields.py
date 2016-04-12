"""
FIRE-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Britton Smith.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.species_fields import \
    add_species_field_by_density
from yt.frontends.gadget.fields import \
    GadgetFieldInfo

class FIREFieldInfo(GadgetFieldInfo):
    def setup_gas_particle_fields(self, ptype):
        super(FIREFieldInfo, self).setup_gas_particle_fields(ptype)
        self.alias((ptype, "temperature"), (ptype, "Temperature"))

        def _h_density(field, data):
            x_H = 1.0 - data[(ptype, "He_fraction")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              data[(ptype, "NeutralHydrogenAbundance")]

        self.add_field(
            (ptype, "H_density"),
            function=_h_density,
            particle_type=True,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H", particle_type=True)
        for suffix in ["density", "fraction", "mass", "number_density"]:
            self.alias((ptype, "H_%s" % suffix), (ptype, "H_p0_%s" % suffix))

        def _h_p1_density(field, data):
            x_H = 1.0 - data[(ptype, "He_fraction")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              (1.0 - data[(ptype, "NeutralHydrogenAbundance")])

        self.add_field(
            (ptype, "H_p1_density"),
            function=_h_p1_density,
            particle_type=True,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H_p1", particle_type=True)
