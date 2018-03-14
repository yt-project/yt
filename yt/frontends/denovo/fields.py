"""
Denovo-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.
#
flux_units = "neutrons / code_length**2"
angular_flux_units = "neutrons / code_length**2 / steradian "
density_units = "code_mass / code_length**3"
source_units = "neutrons"


class DenovoFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("source", (source_units, ["source_strength"], None)),
        ("flux", (flux_units, ["scalar_flux"], r"\phi")),
        ("uncflux", (flux_units, ["uncollided_flux"], r"\phi_{uncollided}")),
        ("angular_flux", (angular_flux_units, [], r"\Psi")),
        ("ww_lower", ("", ["ww_lower_bound"], None)),
    )

    known_particle_fields = (
        # At this time Denovo does not output particle fields, so there are no
        # particle fields defined in this file.
    )

    def __init__(self, ds, field_list):
        super(DenovoFieldInfo, self).__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        params = self.ds.parameters
        self.setup_flux_field()
        pass

    def setup_flux_field(self):
        # Right now I feel that it would be more intuitive to a new user to
        # have flux fields loaded in a descriptive field, however, I might
        # change this if we load in materials information, or ww information,
        # which probably shouldn't be added to this function.
        unit_system = self.ds.unit_system
        params = self.ds.parameters

        super(DenovoFieldInfo, self).setup_flux_fields(fluxtype)

    def setup_particle_fields(self, ptype):
        # Becuase there are no particle types in Denovo, this is empty for now.
        super(DenovoFieldInfo, self).setup_particle_fields(ptype)

