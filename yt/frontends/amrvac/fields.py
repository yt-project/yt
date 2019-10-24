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
from yt.units import dimensions

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.
class AMRVACFieldInfo(FieldInfoContainer):
    code_density = "code_mass / code_length**3"
    code_momentum = "code_mass / code_length**2 / code_time"
    code_energy = "code_mass * code_length**2 / code_time**2"

    # format: (native(?) field, (units, [aliases], display_name))
    # note: aliases will correspond to "gas" typed fields, whereas the native ones are "amrvac" typed
    known_other_fields = (
        ("rho", (code_density, ["density"], r"$\rho$")),
        ("m1", (code_momentum, ["momentum_1"], r"$m_1$")),
        ("m2", (code_momentum, ["momentum_2"], r"$m_2$")),
        ("m3", (code_momentum, ["momentum_3"], r"$m_3$")),
        ("e", (code_energy, ["energy"], r"$e$")),
        ("b1", ("code_magnetic", ["magnetic_1"], r"$B_1$")),
        ("b2", ("code_magnetic", ["magnetic_2"], r"$B_2$")),
        ("b3", ("code_magnetic", ["magnetic_3"], r"$B_3$"))
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        def _v1(field, data):
            return data["gas", "momentum_1"] / data["gas", "density"]
        def _v2(field, data):
            return data["gas", "momentum_2"] / data["gas", "density"]
        def _v3(field, data):
            return data["gas", "momentum_3"] / data["gas", "density"]

        for idir, func in zip("123", (_v1, _v2, _v3)):
            if not ("amrvac", "m%s" % idir) in self.field_list:
                break
            self.add_field(("gas", "velocity_%s" % idir), function=func,
                            display_name=r"$v_%s$" % idir,
                            units="auto", dimensions=dimensions.velocity, sampling_type="cell")

        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])
