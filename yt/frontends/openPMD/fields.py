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


def _kinetic_energy(field, data):
    """
    Function to calculate new fields out of other fields
    if a error occurs the field is not calculated
    """
    #mylog.info("oPMD - fields - _kinetic_energy")
    # calculate kinetic energy out of momentum
    c = 2.997e8  # velocity of light
    return (data["particle_momentum_x"]**2 + data["particle_momentum_y"]**2 + data["particle_momentum_z"]**2) * c


def setup_velocity(self, ptype):
    mylog.info("oPMD - fields - setup_momentum_to_velocity")
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
        mylog.info("oPMD - fields - setup_momentum_to_velocity - particle_velocity_{}".format(ax))
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
            Bfieldx = data["B_x"]
            Bfieldy = data["B_y"]
            Bfieldz = data["B_z"]

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
                       units="T*V/m")  # N/(m*s)


def setup_magnetic(self):
    def _get_mag(axis):
        def magnetic(field, data):
            return data["B_{}".format(axis)]

        return magnetic

    for ax in 'xyz':
        self.add_field(("openPMD", "magnetic_field_%s" % ax),
                       function=_get_mag(ax),
                       units="T")


def setup_electrical(self):
    def _get_el(axis):
        def electrical(field, data):
            return data["E_{}".format(axis)]

        return electrical

    for ax in 'xyz':
        self.add_field(("openPMD", "electrical_field_%s" % ax),
                       function=_get_el(ax),
                       units="kg*m/(A*s**3)")


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

    def __init__(self, ds, field_list):
        super(openPMDFieldInfo, self).__init__(ds, field_list)

        other_fields = ()
        f = ds._handle
        fields = f[ds.basePath + f.attrs["meshesPath"]]
        for i in fields.keys():
            field = fields.get(i)
            for j in field.attrs["axisLabels"]:
                parsed = parse_unitDimension(np.asarray(field.attrs["unitDimension"], dtype='int'))
                other_fields += ((str("_".join([i.replace("_","-"),j])), (yt.YTQuantity(1, parsed).units, [], None)),)
        self.known_other_fields = other_fields
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


    def setup_fluid_fields(self):
        """
        Here you can create functions to calculate out of existing fields the
        values of new fields e.g. calculate out of E-field and B-field the
        Poynting vector
        """
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field and self.add_field .

        setup_poynting_vector(self)
        setup_magnetic(self)
        setup_electrical(self)
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        for ax in "xyz":
            mylog.info("setup_magnetic_field_aliases(self, openPMD, B_{}".format(ax))
            setup_magnetic_field_aliases(self, "openPMD", ["magnetic_field_{}".format(ax)])

    def setup_particle_fields(self, ptype):
        """
        This will get called for every particle type.

        TODO You have to call the function of the parent class to load particles
        """
        # self.add_field(
        #     (ptype,
        #      "particle_kinetic_energy"),
        #     function=_kinetic_energy,
        #     units="dimensionless")
        setup_velocity(self, ptype)
        setup_relative_positions(self, ptype)
        super(openPMDFieldInfo, self).setup_particle_fields(ptype)
