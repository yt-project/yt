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

import numpy as np
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.units import dimensions
from yt import mylog

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.

direction_aliases = {
    "cartesian": ("x", "y", "z"),
    "polar": ("r", "theta", "z"),
    "cylindrical": ("r", "z", "theta"),
    "spherical": ("r", "theta", "phi")
}

def _velocity(idir, idust=None):
    """Generate a velocity function v_idir = m_idir / rho

    If idust is not used, the "gas" velocity function is returned.

    Parameters
    ----------
    idir : int
        the direction index (1, 2 or 3)
    idust : int, optional
        the dust species index (AMRVAC). idust >= 1
    """
    if idust is None:
        dust_label = ""
    else:
        dust_label = "dust%d_" % idust
    def velocity_idir(field, data):
        return data["gas", "%smoment_%d" % (dust_label, idir)] / data["gas", "%sdensity" % dust_label]
    return velocity_idir

code_density = "code_mass / code_length**3"
code_moment = "code_mass / code_length**2 / code_time"
code_pressure = "code_mass / code_length / code_time**2"

# for now, define a finite family of dust fields (up to 100 species, should be enough)
MAXN_DUST_SPECIES = 100
known_dust_fields = [("rhod%d" % idust, (code_density, ["dust%d_density" % idust], None))
                     for idust in range(1, MAXN_DUST_SPECIES+1)]
for idir in range(1, 4):
    known_dust_fields += [("m%dd%d" % (idir, idust), (code_moment, ["dust%d_moment_%d" % (idust, idir)], None))
                          for idust in range(1, MAXN_DUST_SPECIES+1)]

class AMRVACFieldInfo(FieldInfoContainer):

    # format: (native(?) field, (units, [aliases], display_name))
    # note: aliases will correspond to "gas" typed fields, whereas the native ones are "amrvac" typed
    known_other_fields = tuple(list((
        ("rho", (code_density, ["density"], None)),
        ("m1", (code_moment, ["moment_1"], None)),
        ("m2", (code_moment, ["moment_2"], None)),
        ("m3", (code_moment, ["moment_3"], None)),
        ("e", (code_pressure, ["energy_density"], None)),
        ("b1", ("code_magnetic", ["magnetic_1"], None)),
        ("b2", ("code_magnetic", ["magnetic_2"], None)),
        ("b3", ("code_magnetic", ["magnetic_3"], None))
        )) + known_dust_fields
        # in python3, there is no need for this tuple+list conversion, it suffices to write
        # known_other_fields = (..., *known_dust_fields)
    )

    known_particle_fields = ()

    def _create_velocity_fields(self, idust=None):
        if idust is None:
            dust_flag = dust_label = ""
        else:
            dust_flag = "d%d" % idust
            dust_label = "dust%d_" % idust

        unit_system = self.ds.unit_system
        for i_, alias in enumerate(direction_aliases[self.ds.geometry]):
            idir = i_ + 1 # direction index starts at 1 in AMRVAC
            if not ("amrvac", "m%d%s" % (idir, dust_flag)) in self.field_list:
                break
            self.add_field(("gas", "%svelocity_%s" % (dust_label, alias)), function=_velocity(idir, idust),
                            units=unit_system['velocity'], dimensions=dimensions.velocity, sampling_type="cell")
            self.alias(("gas", "%svelocity_%d" % (dust_label, idir)), ("gas", "%svelocity_%s" % (dust_label, alias)),
                        units=unit_system["velocity"])
            self.alias(("gas", "%smoment_%s" % (dust_label, alias)), ("gas", "%smoment_%d" % (dust_label, idir)),
                        units=unit_system["density"]*unit_system["velocity"])

    def setup_fluid_fields(self):

        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])
        self._create_velocity_fields() # gas velocity

        # dust derived fields
        unit_system = self.ds.unit_system
        idust = 1
        while ("amrvac", "rhod%d" % idust) in self.field_list:
            if idust > MAXN_DUST_SPECIES:
                mylog.error("Only the first %d dust species are currently read by yt. " \
                            "If you read this, please consider issuing a ticket. " % MAXN_DUST_SPECIES)
                break
            self._create_velocity_fields(idust)
            idust += 1
        n_dust_found = idust - 1

        def _total_dust_density(field, data):
            tot = np.zeros_like(data["density"])
            for idust in range(1, n_dust_found+1):
                tot += data["dust%d_density" % idust]
            return tot

        self.add_field(("gas", "total_dust_density"), function=_total_dust_density,
                       dimensions=dimensions.density,
                       units=unit_system["density"], sampling_type="cell")

        def dust_to_gas_ratio(field, data):
            return data["total_dust_density"] / data["density"]

        self.add_field(("gas", "dust_to_gas_ratio"), function=dust_to_gas_ratio,
                        dimensions=dimensions.dimensionless,
                        sampling_type="cell")

