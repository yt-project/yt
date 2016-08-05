"""
openPMD-specific fields



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from yt.funcs import mylog
from yt.utilities.physical_constants import speed_of_light
from yt.fields.field_info_container import FieldInfoContainer
from yt.units.yt_array import YTQuantity
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.frontends.open_pmd.misc import parse_unit_dimension

import numpy as np


def setup_kinetic_energy(self, ptype):
    def _kin_en(field, data):
        p2 = (data[ptype, "particle_momentum_x"] ** 2 +
              data[ptype, "particle_momentum_y"] ** 2 +
              data[ptype, "particle_momentum_z"] ** 2)
        mass = data[ptype, "particle_mass"] * data[ptype, "particle_weighting"]
        return speed_of_light * np.sqrt(p2 + mass ** 2 * speed_of_light ** 2) - mass * speed_of_light ** 2

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
                                  (mass * weighting) ** 2 +
                                  (momentum ** 2) / (c ** 2)
                              ) ** 0.5

        return velocity

    for ax in "xyz":
        self.add_field((ptype, "particle_velocity_%s" % ax),
                       function=_get_vel(ax),
                       units="m/s",
                       particle_type=True)


def setup_poynting_vector(self):
    def _get_poyn(axis):
        def poynting(field, data):
            u = 79577.4715459  # = 1/magnetic permeability
            if axis in "x":
                return u * (data["E_y"] * data["magnetic_field_z"] - data["E_z"] * data["magnetic_field_y"])
            elif axis in "y":
                return u * (data["E_z"] * data["magnetic_field_x"] - data["E_x"] * data["magnetic_field_z"])
            elif axis in "z":
                return u * (data["E_x"] * data["magnetic_field_y"] - data["E_y"] * data["magnetic_field_x"])

        return poynting

    for ax in "xyz":
        self.add_field(("openPMD", "poynting_vector_%s" % ax),
                       function=_get_poyn(ax),
                       units="T*V/m")


def setup_absolute_positions(self, ptype):
    def _abs_pos(axis):
        def ap(field, data):
            return np.add(data[ptype, "particle_positionCoarse_{}".format(axis)],
                          data[ptype, "particle_positionOffset_{}".format(axis)])

        return ap

    for ax in "xyz":
        self.add_field((ptype, "particle_position_%s" % ax),
                       function=_abs_pos(ax),
                       units="m",
                       particle_type=True)


class OpenPMDFieldInfo(FieldInfoContainer):
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
        bp = ds.base_path
        mp = ds.meshes_path
        pp = ds.particles_path
        fields = f[bp + mp]

        for fname in fields.keys():
            field = fields.get(fname)
            if "dataset" in str(field).split(" ")[1]:
                # We have a dataset, don't consider axes. This appears to be a vector field of single dimensionality
                ytname = str("_".join([fname.replace("_", "-")]))
                if ds._nonstandard:
                    parsed = ""
                else:
                    parsed = parse_unit_dimension(np.asarray(field.attrs["unitDimension"], dtype="int"))
                unit = str(YTQuantity(1, parsed).units)
                aliases = []
                # Save a list of magnetic fields for aliasing later on
                # We can not reasonably infer field type by name in openPMD
                if "T" in unit or "kg/(A*s**2)" in unit:
                    self._mag_fields.append(ytname)
                self.known_other_fields += ((ytname, (unit, aliases, None)),)
            else:
                if ds._nonstandard:
                    axes = "xyz"  # naively assume all fields in non-standard files are 3D
                else:
                    axes = field.attrs["axisLabels"]
                for axis in axes:
                    ytname = str("_".join([fname.replace("_", "-"), axis]))
                    if ds._nonstandard:
                        parsed = ""
                    else:
                        parsed = parse_unit_dimension(np.asarray(field.attrs["unitDimension"], dtype="int"))
                    unit = str(YTQuantity(1, parsed).units)
                    aliases = []
                    # Save a list of magnetic fields for aliasing later on
                    # We can not reasonably infer field type by name in openPMD
                    if "T" in unit or "kg/(A*s**2)" in unit:
                        self._mag_fields.append(ytname)
                    self.known_other_fields += ((ytname, (unit, aliases, None)),)
        for i in self.known_other_fields:
            mylog.debug("oPMD - fields - known_other_fields - {}".format(i))

        particle_fields = ()
        particles = f[bp + pp]
        for species in particles.keys():
            for attrib in particles.get(species).keys():
                if "weighting" in attrib:
                    particle_fields += (("particle_weighting", ("", [], None)),)
                    continue
                try:
                    if ds._nonstandard:
                        if "globalCellIdx" in attrib or "position" in attrib:
                            parsed = "m"  # Required for spatial selection of particles
                        else:
                            parsed = ""
                    else:
                        parsed = parse_unit_dimension(
                            np.asarray(particles.get(species).get(attrib).attrs["unitDimension"], dtype="int"))
                    unit = str(YTQuantity(1, parsed).units)
                    name = ["particle", attrib]
                    ytattrib = attrib
                    # Symbolically rename position to preserve yt's interpretation of the pfield
                    # particle_position is later derived in setup_absolute_positions
                    if ytattrib in "position":
                        ytattrib = "positionCoarse"
                    for axis in particles.get(species).get(attrib).keys():
                        aliases = []
                        if axis in "rxyz":
                            name = ["particle", ytattrib, axis]
                        ytname = str("_".join(name))
                        if ds._nonstandard and "globalCellIdx" in ytname:
                            aliases.append(ytname.replace("globalCellIdx", "positionOffset"))
                        particle_fields += ((ytname, (unit, aliases, None)),)
                except:
                    mylog.info("{}_{} does not seem to have unitDimension".format(species, attrib))
        self.known_particle_fields = particle_fields
        for i in self.known_particle_fields:
            mylog.debug("oPMD - fields - known_particle_fields - {}".format(i))
        super(OpenPMDFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        """
        Here you can create functions to calculate out of existing fields the
        values of new fields e.g. calculate out of E-field and B-field the
        Poynting vector
        """
        # Set up aliases first so the setup for poynting can use them
        if len(self._mag_fields) > 0:
            setup_magnetic_field_aliases(self, "openPMD", self._mag_fields)
            setup_poynting_vector(self)

    def setup_particle_fields(self, ptype):
        """
        This will get called for every particle type.
        """
        setup_absolute_positions(self, ptype)
        setup_kinetic_energy(self, ptype)
        setup_velocity(self, ptype)
        super(OpenPMDFieldInfo, self).setup_particle_fields(ptype)
