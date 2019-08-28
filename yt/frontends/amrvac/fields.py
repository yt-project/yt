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
b_units   = "code_magnetic"

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class AMRVACFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("rho", (rho_units, ["rho", "to", "density"], r"$\rho$")),
        ("m1", (mom_units, ["momentum_1"], r"$m_1$")),
        ("m2", (mom_units, ["momentum_2"], r"$m_2$")),
        ("m3", (mom_units, ["momentum_3"], r"$m_3$")),
        ("e", ("J", ["energy"], r"$e$")),
        ("b1", (b_units, [], r"$B_x$")),
        ("b2", (b_units, [], r"$B_y$")),
        ("b3", (b_units, [], r"$B_z$"))
    )

    known_particle_fields = ()

    def __init__(self, ds, field_list):
        super(AMRVACFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases

        # Add pressure
        # def _get_ekin(data):
        #     ekin = 0.5 * data["amrvac", "m1"]**2 / data["amrvac", "rho"]
        #     if self.ds.dimensionality > 1:
        #         ekin = ekin + 0.5 * data["amrvac", "m2"] ** 2 / data["amrvac", "rho"]
        #     if self.ds.dimensionality > 2:
        #         ekin = ekin + 0.5 * data["amrvac", "m3"] ** 2 / data["amrvac", "rho"]
        #     return ekin
        # def _pressure(field, data):
        #     return (self.ds.gamma - 1) * (data["amrvac", "e"] - _get_ekin(data))

        # self.add_field(("amrvac", "pressure"), sampling_type="cell",
        #                function=_pressure, units="code_pressure")
        # self.alias(("amrvac", "rho"), ("amrvac", "p"), units="code_pressure")
        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])
        # print(self[("amrvac", "rho")])
        # print(self[("amrvac", "pressure")])


    def setup_particle_fields(self, ptype):
        super(AMRVACFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
