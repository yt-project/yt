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
import yt.utilities.physical_constants
from yt.fields.field_info_container import \
    FieldInfoContainer, \
    particle_deposition_functions, \
    particle_vector_functions, \
    standard_particle_fields

def _kinetic_energy(field, data):
    """
    Function to calculate new fields out of other fields
    if a error occurs the field is not calculated
    """
    mylog.info("oPMD - fields - _kinetic_energy")
    # calculate kinetic energy out of momentum
    c = 2.997e8  # velocity of light
    ke = (data["particle_momentum_x"] ** 2
          + data["particle_momentum_y"] ** 2
          + data["particle_momentum_z"] ** 2) * c
    return ke


def setup_momentum_to_velocity(self, ptype):
    """
    TODO This function does no work !!
    """
    mylog.info("oPMD - fields - setup_momentum_to_velocity")
    def _get_vel(axis):
        def velocity(field, data):
            c = 2.997e8
            moment = data[ptype, "particle_momentum_%" % axis]
            return moment / ((data[ptype, "particle_mass"] * data[ptype, "particle_weighting"]) ** 2 + (moment ** 2) / (c ** 2)) ** 0.5
        return velocity

    for ax in 'xyz':
        self.add_field((ptype, "particle_velocity_%s" % ax),
                       function=_get_vel(ax),
                       units="m/s")


def setup_poynting_vector(self):
    mylog.info("oPMD - fields - setup_poynting_vector")
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


class openPMDFieldInfo(FieldInfoContainer):
    """
    We need to specify which fields we might have in our dataset.  The field info
    container subclass here will define which fields it knows about.  There are
    optionally methods on it that get called which can be subclassed.

    !! TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    This class defines, which fields and particle fields could be in the HDF5-file
    The field names have to match the names in "openPMDHierarchy" in data_structures.py
    This also defines the units of the fields
    """
    # TODO Generate this when parsing a file
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("B_x", ("T", [], None)),
        ("B_y", ("T", [], None)),
        ("B_z", ("T", [], None)),
        ("E_x", ("V/m", [], None)),
        ("E_y", ("V/m", [], None)),
        ("E_z", ("V/m", [], None)),
        ("J_x", ("A/m**2", [], None)),
        ("J_y", ("A/m**2", [], None)),
        ("J_z", ("A/m**2", [], None)),
        ("rho", ("kg*s/m**3.0", [], None)),

    )

    # TODO Generate this when parsing a file
    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("particle_charge", ("A*s", [], None)),
        ("particle_mass", ("kg", [], None)),
        ("particle_momentum_x", ("kg*m/s", [], None)),
        ("particle_momentum_y", ("kg*m/s", [], None)),
        ("particle_momentum_z", ("kg*m/s", [], None)),
        ("particle_position_x", ("m", [], None)),
        ("particle_position_y", ("m", [], None)),
        ("particle_position_z", ("m", [], None)),
        ("particle_positionOffset_x", ("m", [], None)),
        ("particle_positionOffset_y", ("m", [], None)),
        ("particle_positionOffset_z", ("m", [], None)),
        ("particle_weighting", ("", [], None)),
        ("particle_kinetic_energy", ("dimensionless", [], None)),

    )

    def __init__(self, ds, field_list):
        super(openPMDFieldInfo, self).__init__(ds, field_list)
        # If you want, you can check field_list

    def setup_fluid_fields(self):
        """
        Here you can create functions to calculate out of existing fields the
        values of new fields e.g. calculate out of E-field and B-field the
        Poynting vector
        """
        mylog.info("oPMD - fields - setup_fluid_fields")
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field and self.add_field .

        setup_poynting_vector(self)

    def setup_particle_fields(self, ptype):
        """
        This will get called for every particle type.

        TODO You have to call the function of the parent class to load particles
        """

        mylog.info("oPMD - fields - setup_particle_fields(%s)",ptype)
        self.add_field(
            (ptype,
             "particle_kinetic_energy"),
            function=_kinetic_energy,
            units="dimensionless")
        setup_momentum_to_velocity(self, ptype)

        particle_deposition_functions(ptype, "particle_position",
                                      "particle_weighting", self)
        standard_particle_fields(self, ptype)

        # TODO Has to be called to load particles
        super(openPMDFieldInfo, self).setup_particle_fields(ptype)
