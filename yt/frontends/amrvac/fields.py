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

direction_aliases = {
    "cartesian": ("x", "y", "z"),
    "polar": ("r", "theta", "z"),
    "cylindrical": ("r", "z", "theta"),
    "spherical": ("r", "theta", "phi")
}

class AMRVACFieldInfo(FieldInfoContainer):
    code_density = "code_mass / code_length**3"
    code_moment = "code_mass / code_length**2 / code_time"
    code_pressure = "code_mass / code_length / code_time**2"

    # format: (native(?) field, (units, [aliases], display_name))
    # note: aliases will correspond to "gas" typed fields, whereas the native ones are "amrvac" typed

    # for now, define a finite family of dust fields (up to 100 species, should be enough)
    dust_fields = [("rhod%d" % idust, ("code_mass / code_length**3", ["dust%d_density" % idust], None))
                   for idust in range(1, 100)]
    known_other_fields = (
        ("rho", (code_density, ["density"], None)),
        ("m1", (code_moment, ["moment_1"], None)),
        ("m2", (code_moment, ["moment_2"], None)),
        ("m3", (code_moment, ["moment_3"], None)),
        ("e", (code_pressure, ["energy_density"], None)),
        ("b1", ("code_magnetic", ["magnetic_1"], None)),
        ("b2", ("code_magnetic", ["magnetic_2"], None)),
        ("b3", ("code_magnetic", ["magnetic_3"], None)),
        *dust_fields
    )


    known_particle_fields = ()

    def setup_fluid_fields(self):

        unit_system = self.ds.unit_system
        aliases = direction_aliases[self.ds.geometry]

        def _velocity(idir):
            """Generate a velocity function v_idir = m_idir / rho

            Parameters
            ----------
            idir : int
                the direction index (1, 2 or 3)
            """
            def velocity_idir(field, data):
                return data["gas", "moment_%d" % idir] / data["gas", "density"]
            return velocity_idir

        for i, alias in enumerate(aliases):
            idir = i + 1 # direction index starts at 1 in AMRVAC
            if not ("amrvac", "m%s" % idir) in self.field_list:
                break
            self.add_field(("gas", "velocity_%s" % alias), function=_velocity(idir),
                            units=unit_system['velocity'], dimensions=dimensions.velocity, sampling_type="cell")
            self.alias(("gas", "velocity_%d" % idir), ("gas", "velocity_%s" % alias),
                        units=unit_system["velocity"])
            self.alias(("gas", "moment_%s" % alias), ("gas", "moment_%d" % idir),
                        units=unit_system["density"]*unit_system["velocity"])

        # dust
        idust = 1
        while ("amrvac", "rhod%d" % idust) in self.field_list:
        #    self.add_field(("gas", "dust%d_density" % idust), function=_density(idust),
        #                    units=unit_system['density'], dimensions=dimensions.density, sampling_type="cell")
            idust += 1


        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])
