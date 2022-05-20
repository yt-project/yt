from typing import List

import numpy as np

from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.frontends.open_pmd.misc import is_const_component, parse_unit_dimension
from yt.units.yt_array import YTQuantity
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.physical_constants import mu_0, speed_of_light


def setup_poynting_vector(self):
    def _get_poyn(axis):
        def poynting(field, data):
            u = mu_0**-1
            if axis in "x":
                return u * (
                    data[("openPMD", "E_y")] * data[("gas", "magnetic_field_z")]
                    - data[("openPMD", "E_z")] * data[("gas", "magnetic_field_y")]
                )
            elif axis in "y":
                return u * (
                    data[("openPMD", "E_z")] * data[("gas", "magnetic_field_x")]
                    - data[("openPMD", "E_x")] * data[("gas", "magnetic_field_z")]
                )
            elif axis in "z":
                return u * (
                    data[("openPMD", "E_x")] * data[("gas", "magnetic_field_y")]
                    - data[("openPMD", "E_y")] * data[("gas", "magnetic_field_x")]
                )

        return poynting

    for ax in "xyz":
        self.add_field(
            ("openPMD", f"poynting_vector_{ax}"),
            sampling_type="cell",
            function=_get_poyn(ax),
            units="W/m**2",
        )


def setup_kinetic_energy(self, ptype):
    def _kin_en(field, data):
        p2 = (
            data[ptype, "particle_momentum_x"] ** 2
            + data[ptype, "particle_momentum_y"] ** 2
            + data[ptype, "particle_momentum_z"] ** 2
        )
        mass = data[ptype, "particle_mass"] * data[ptype, "particle_weighting"]
        return (
            speed_of_light * np.sqrt(p2 + mass**2 * speed_of_light**2)
            - mass * speed_of_light**2
        )

    self.add_field(
        (ptype, "particle_kinetic_energy"),
        sampling_type="particle",
        function=_kin_en,
        units="kg*m**2/s**2",
    )


def setup_velocity(self, ptype):
    def _get_vel(axis):
        def velocity(field, data):
            c = speed_of_light
            momentum = data[ptype, f"particle_momentum_{axis}"]
            mass = data[ptype, "particle_mass"]
            weighting = data[ptype, "particle_weighting"]
            return momentum / np.sqrt(
                (mass * weighting) ** 2 + (momentum**2) / (c**2)
            )

        return velocity

    for ax in "xyz":
        self.add_field(
            (ptype, f"particle_velocity_{ax}"),
            sampling_type="particle",
            function=_get_vel(ax),
            units="m/s",
        )


def setup_absolute_positions(self, ptype):
    def _abs_pos(axis):
        def ap(field, data):
            return np.add(
                data[ptype, f"particle_positionCoarse_{axis}"],
                data[ptype, f"particle_positionOffset_{axis}"],
            )

        return ap

    for ax in "xyz":
        self.add_field(
            (ptype, f"particle_position_{ax}"),
            sampling_type="particle",
            function=_abs_pos(ax),
            units="m",
        )


class OpenPMDFieldInfo(FieldInfoContainer):
    r"""Specifies which fields from the dataset yt should know about.

    ``self.known_other_fields`` and ``self.known_particle_fields`` must be populated.
    Entries for both of these lists must be tuples of the form ("name", ("units",
    ["fields", "to", "alias"], "display_name")) These fields will be represented and
    handled in yt in the way you define them here. The fields defined in both
    ``self.known_other_fields`` and ``self.known_particle_fields`` will only be added to
    a dataset (with units, aliases, etc), if they match any entry in the
    ``OpenPMDHierarchy``'s ``self.field_list``.

    Notes
    -----

    Contrary to many other frontends, we dynamically obtain the known fields from the
    simulation output. The openPMD markup is extremely flexible - names, dimensions and
    the number of individual datasets can (and very likely will) vary.

    openPMD states that record names and their components are only allowed to contain
    * characters a-Z,
    * the numbers 0-9
    * and the underscore _
    * (equivalently, the regex \w).
    Since yt widely uses the underscore in field names, openPMD's underscores (_) are
    replaced by hyphen (-).

    Derived fields will automatically be set up, if names and units of your known
    on-disk (or manually derived) fields match the ones in [1].

    References
    ----------
    * http://yt-project.org/docs/dev/analyzing/fields.html
    * http://yt-project.org/docs/dev/developing/creating_frontend.html#data-meaning-structures
    * https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md
    * [1] http://yt-project.org/docs/dev/reference/field_list.html#universal-fields
    """

    _mag_fields: List[str] = []

    def __init__(self, ds, field_list):
        f = ds._handle
        bp = ds.base_path
        mp = ds.meshes_path
        pp = ds.particles_path

        try:
            fields = f[bp + mp]
            for fname in fields.keys():
                field = fields[fname]
                if isinstance(field, h5py.Dataset) or is_const_component(field):
                    # Don't consider axes.
                    # This appears to be a vector field of single dimensionality
                    ytname = str("_".join([fname.replace("_", "-")]))
                    parsed = parse_unit_dimension(
                        np.asarray(field.attrs["unitDimension"], dtype="int64")
                    )
                    unit = str(YTQuantity(1, parsed).units)
                    aliases = []
                    # Save a list of magnetic fields for aliasing later on
                    # We can not reasonably infer field type/unit by name in openPMD
                    if unit == "T" or unit == "kg/(A*s**2)":
                        self._mag_fields.append(ytname)
                    self.known_other_fields += ((ytname, (unit, aliases, None)),)
                else:
                    for axis in field.keys():
                        ytname = str("_".join([fname.replace("_", "-"), axis]))
                        parsed = parse_unit_dimension(
                            np.asarray(field.attrs["unitDimension"], dtype="int64")
                        )
                        unit = str(YTQuantity(1, parsed).units)
                        aliases = []
                        # Save a list of magnetic fields for aliasing later on
                        # We can not reasonably infer field type by name in openPMD
                        if unit == "T" or unit == "kg/(A*s**2)":
                            self._mag_fields.append(ytname)
                        self.known_other_fields += ((ytname, (unit, aliases, None)),)
            for i in self.known_other_fields:
                mylog.debug("open_pmd - known_other_fields - %s", i)
        except (KeyError, TypeError, AttributeError):
            pass

        try:
            particles = f[bp + pp]
            for pname in particles.keys():
                species = particles[pname]
                for recname in species.keys():
                    try:
                        record = species[recname]
                        parsed = parse_unit_dimension(record.attrs["unitDimension"])
                        unit = str(YTQuantity(1, parsed).units)
                        ytattrib = str(recname).replace("_", "-")
                        if ytattrib == "position":
                            # Symbolically rename position to preserve yt's
                            # interpretation of the pfield particle_position is later
                            # derived in setup_absolute_positions in the way yt expects
                            ytattrib = "positionCoarse"
                        if isinstance(record, h5py.Dataset) or is_const_component(
                            record
                        ):
                            name = ["particle", ytattrib]
                            self.known_particle_fields += (
                                (str("_".join(name)), (unit, [], None)),
                            )
                        else:
                            for axis in record.keys():
                                aliases = []
                                name = ["particle", ytattrib, axis]
                                ytname = str("_".join(name))
                                self.known_particle_fields += (
                                    (ytname, (unit, aliases, None)),
                                )
                    except (KeyError):
                        if recname != "particlePatches":
                            mylog.info(
                                "open_pmd - %s_%s does not seem to have "
                                "unitDimension",
                                pname,
                                recname,
                            )
            for i in self.known_particle_fields:
                mylog.debug("open_pmd - known_particle_fields - %s", i)
        except (KeyError, TypeError, AttributeError):
            pass

        super().__init__(ds, field_list)

    def setup_fluid_fields(self):
        """Defines which derived mesh fields to create.

        If a field can not be calculated, it will simply be skipped.
        """
        # Set up aliases first so the setup for poynting can use them
        if len(self._mag_fields) > 0:
            setup_magnetic_field_aliases(self, "openPMD", self._mag_fields)
            setup_poynting_vector(self)

    def setup_particle_fields(self, ptype):
        """Defines which derived particle fields to create.

        This will be called for every entry in
        `OpenPMDDataset``'s ``self.particle_types``.
        If a field can not be calculated, it will simply be skipped.
        """
        setup_absolute_positions(self, ptype)
        setup_kinetic_energy(self, ptype)
        setup_velocity(self, ptype)
        super().setup_particle_fields(ptype)
