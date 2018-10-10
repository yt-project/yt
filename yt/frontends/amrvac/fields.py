"""
AMRVAC-specific fields



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

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_length**2*code_time)"

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class AMRVACFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("density", (rho_units, ["density"], r"$\rho$")),
        ("m1", (mom_units, ["momentum_1"], r"$m_1$")),
        ("m2", (mom_units, ["momentum_2"], r"$m_2$")),
        ("m3", (mom_units, ["momentum_3"], r"$m_3$")),
        ("energy", ("J", ["energy"], r"$e$")),
        #devnote : missing a way to handle an arbitrary number of dust fluids here
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def __init__(self, ds, field_list):
        super(AMRVACFieldInfo, self).__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        pass

    def setup_particle_fields(self, ptype):
        super(AMRVACFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
