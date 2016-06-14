"""
openPMD-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog
from yt.utilities.physical_constants import speed_of_light
from yt.fields.field_info_container import \
    FieldInfoContainer
from .misc import parse_unitDimension
import yt


def setup_kinetic_energy(self, ptype):
    def _kin_en(field, data):
        # Calculation seems to be wrong:
        # YTFieldUnitError
        # The field function associated with the field '('io', 'particle_kinetic_energy')'
        # returned data with units 'cm*kg**2*m**2/s**3' but was defined with units 'kg*m**2/s**2'.
        return (
                   data[ptype, "particle_momentum_x"]**2 +
                   data[ptype, "particle_momentum_y"]**2 +
                   data[ptype, "particle_momentum_z"]**2
               ) * speed_of_light

    self.add_field((ptype, "particle_kinetic_energy"),
                   function=_kin_en,
                   units="kg*m**2/s**2",
                   particle_type=True)


def setup_velocity(self, ptype):
    def _get_vel(axis):
        def velocity(field, data):
            c = speed_of_light
            momentum = data[ptype, "particle_momentum_{}".format(axis)]
            mass = data[ptype, "particle_mass"]
            weighting = data[ptype, "particle_weighting"]
            return momentum / (
                                  (mass * weighting)**2 +
                                  (momentum**2) / (c**2)
                              ) ** 0.5

        return velocity

    for ax in 'xyz':
        self.add_field((ptype, "particle_velocity_%s" % ax),
                       function=_get_vel(ax),
                       units="m/s",
                       particle_type=True)


def setup_poynting_vector(self):
    def _get_poyn(axis):
        def poynting(field, data):
            Efieldx = data["E_x"]
            Efieldy = data["E_y"]
            Efieldz = data["E_z"]
            Bfieldx = data["magnetic_field_x"]
            Bfieldy = data["magnetic_field_y"]
            Bfieldz = data["magnetic_field_z"]

            u = 79577.4715459  # = 1/magnetic permeability

            if(axis == 'x'):
                return u * (Efieldy * Bfieldz - Efieldz * Bfieldy)
            elif(axis == 'y'):
                return u * (Efieldz * Bfieldx - Efieldx * Bfieldz)
            elif(axis == 'z'):
                return u * (Efieldx * Bfieldy - Efieldy * Bfieldx)

        return poynting

    for ax in 'xyz':
        self.add_field(("openPMD", "poynting_vector_%s" % ax),
                       function=_get_poyn(ax),
                       units="T*V/m")


def setup_relative_positions(self, ptype):
    def _rel_pos(axis):
        def rp(field, data):
            pos = data[ptype, "particle_position_{}".format(axis)]
            off = data[ptype, "particle_positionOffset_{}".format(axis)]
            return pos + off

        return rp

    for ax in 'xyz':
        self.add_field((ptype, "particle_position_relative_%s" % ax),
                       function=_rel_pos(ax),
                       units="m",
                       particle_type=True)


class openPMDFieldInfo(FieldInfoContainer):
    """
    We need to specify which fields we might have in our dataset.  The field info
    container subclass here will define which fields it knows about.  There are
    optionally methods on it that get called which can be subclassed.

    This class defines, which fields and particle fields could be in the HDF5-file
    The field names have to match the names in "openPMDHierarchy" in data_structures.py
    This also defines the units of the fields
    """

    _mag_fields = []

    def __init__(self, ds, field_list):
        f = ds._handle
        fields = f[ds.basePath + f.attrs["meshesPath"]]
        for fname in fields.keys():
            field = fields.get(fname)
            for axis in field.attrs["axisLabels"]:
                ytname = str("_".join([fname.replace("_", "-"), axis]))
                parsed = parse_unitDimension(np.asarray(field.attrs["unitDimension"], dtype='int'))
                unit = str(yt.YTQuantity(1, parsed).units)
                aliases = []
                # Save a list of magnetic fields for aliasing later on
                # We can not reasonably infer field type by name in openPMD
                if "T" in unit or "kg/(A*s**2)" in unit:
                    self._mag_fields.append(ytname)
                self.known_other_fields += ((ytname, (unit, aliases, None)), )
        for i in self.known_other_fields:
            mylog.debug("oPMD - fields - known_other_fields - {}".format(i))

        particle_fields = ()
        particles = f[ds.basePath + f.attrs["particlesPath"]]
        for species in particles.keys():
            for attrib in particles.get(species).keys():
                try:
                    if "weighting" in attrib:
                        particle_fields += (("particle_weighting", ("", [], None)),)
                        continue
                    udim = particles.get(species).get(attrib).attrs["unitDimension"]
                    parsed = parse_unitDimension(np.asarray(udim, dtype='int'))
                    for axis in particles.get(species).get(attrib).keys():
                        if axis in "rxyz":
                            particle_fields += (
                                (str("_".join(["particle", attrib, axis])), (yt.YTQuantity(1, parsed).units, [], None)),)
                        else:
                            particle_fields += (
                            (str("_".join(["particle", attrib])), (yt.YTQuantity(1, parsed).units, [], None)),)
                except:
                    mylog.info("{}_{} does not seem to have unitDimension".format(species, attrib))
        self.known_particle_fields = particle_fields
        for i in self.known_particle_fields:
            mylog.debug("oPMD - fields - known_particle_fields - {}".format(i))
        super(openPMDFieldInfo, self).__init__(ds, field_list)


    def setup_fluid_fields(self):
        """
        Here you can create functions to calculate out of existing fields the
        values of new fields e.g. calculate out of E-field and B-field the
        Poynting vector
        """
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        setup_magnetic_field_aliases(self, "openPMD", self._mag_fields)
        # Set up aliases first so the setup for poynting can use them
        setup_poynting_vector(self)

    def setup_particle_fields(self, ptype):
        """
        This will get called for every particle type.

        TODO You have to call the function of the parent class to load particles
        """
        setup_relative_positions(self, ptype)
        #setup_kinetic_energy(self, ptype)
        setup_velocity(self, ptype)
        super(openPMDFieldInfo, self).setup_particle_fields(ptype)
