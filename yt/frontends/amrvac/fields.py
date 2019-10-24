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
from yt.fields.magnetic_field import setup_magnetic_field_aliases


# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.
class AMRVACFieldInfo(FieldInfoContainer):
    code_density = "code_mass / code_length**3"
    code_momentum = "code_density * code_velocity"

    # format: (native(?) field, (units, [aliases], display_name))
    known_other_fields = (
        ("rho", (code_density, ["density"], r"$\rho$")),
        ("m1", (code_momentum, ["momentum_1"], r"$m_1$")),
        ("m2", (code_momentum, ["momentum_2"], r"$m_2$")),
        ("m3", (code_momentum, ["momentum_3"], r"$m_3$")),
        ("e", ("code_energy", ["energy"], r"$e$")),
        ("b1", ("code_margnetic", [], r"$B_1$")),
        ("b2", ("code_margnetic", [], r"$B_2$")),
        ("b3", ("code_margnetic", [], r"$B_3$"))
    )

    known_particle_fields = ()

    def __init__(self, ds, field_list):
        super(AMRVACFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        # add primitive variables as custom fields
        def _velocity1(field, data):
            return data["amrvac", "m1"] / data["amrvac", "rho"]
        def _velocity2(field, data):
            return data["amrvac", "m2"] / data["amrvac", "rho"]
        def _velocity3(field, data):
            return data["amrvac", "m3"] / data["amrvac", "rho"]

        # velocity fields
        self.add_field(("amrvac", "velocity_1"), sampling_type="cell",
                       function=_velocity1, units="")
        self.alias(("amrvac", "v1"), ("amrvac", "velocity_1"), units="") # missing units

        if ("amrvac", "m2") in self.field_list:
            self.add_field(("amrvac", "velocity_2"), sampling_type="cell",
                           function=_velocity2, units="") # missing units
            self.alias(("amrvac", "v2"), ("amrvac", "velocity_2"), units="") # missing units
        if ("amrvac", "m3") in self.field_list:
            self.add_field(("amrvac", "velocity_3"), sampling_type="cell",
                           function=_velocity3, units="") # missing units
            self.alias(("amrvac", "v3"), ("amrvac", "velocity_3"), units="") # missing units

        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])


    def setup_particle_fields(self, ptype):
        super(AMRVACFieldInfo, self).setup_particle_fields(ptype)
