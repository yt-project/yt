from functools import partial

from yt.fields.particle_fields import sph_whitelist_fields
from yt.frontends.gadget.definitions import elem_names_opts
from yt.frontends.sph.fields import SPHFieldInfo
from yt.utilities.periodic_table import periodic_table
from yt.utilities.physical_constants import kb, mp
from yt.utilities.physical_ratios import _primordial_mass_fraction


class GadgetFieldInfo(SPHFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        if ds.gen_hsmls:
            hsml = (("smoothing_length", ("code_length", [], None)),)
            self.known_particle_fields += hsml
        for field in field_list:
            if field[1].startswith("MetalMasses"):
                mm = ((field[1], ("code_mass", [], None)),)
                self.known_particle_fields += mm
        super().__init__(ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):
        # setup some special fields that only make sense for SPH particles

        if (ptype, "FourMetalFractions") in self.ds.field_list:
            self.species_names = self._setup_four_metal_fractions(ptype)
        elif (ptype, "ElevenMetalMasses") in self.ds.field_list:
            self.species_names = self._setup_metal_masses(ptype, "ElevenMetalMasses")
        elif (ptype, "MetalMasses_00") in self.ds.field_list:
            self.species_names = self._setup_metal_masses(ptype, "MetalMasses")
        if len(self.species_names) == 0:
            self.species_names = self._check_whitelist_species_fields(ptype)

        super().setup_particle_fields(ptype, *args, **kwargs)

        if ptype in ("PartType0", "Gas"):
            self.setup_gas_particle_fields(ptype)

    def _setup_four_metal_fractions(self, ptype):
        """
        This function breaks the FourMetalFractions field (if present)
        into its four component metal fraction fields and adds
        corresponding metal density fields which will later get smoothed

        This gets used with the Gadget group0000 format
        as defined in the gadget_field_specs in frontends/gadget/definitions.py
        """
        metal_names = elem_names_opts[4]

        def _fraction(field, data, i: int):
            return data[(ptype, "FourMetalFractions")][:, i]

        def _metal_density(field, data, i: int):
            return data[(ptype, "FourMetalFractions")][:, i] * data[(ptype, "density")]

        for i, metal_name in enumerate(metal_names):
            # add the metal fraction fields
            self.add_field(
                (ptype, metal_name + "_fraction"),
                sampling_type="particle",
                function=partial(_fraction, i=i),
                units="",
            )

            # add the metal density fields
            self.add_field(
                (ptype, metal_name + "_density"),
                sampling_type="particle",
                function=partial(_metal_density, i=i),
                units=self.ds.unit_system["density"],
            )

        return metal_names

    def _make_fraction_functions(self, ptype, fname):
        if fname == "ElevenMetalMasses":

            def _fraction(field, data, i: int):
                return (
                    data[(ptype, "ElevenMetalMasses")][:, i]
                    / data[(ptype, "particle_mass")]
                )

            def _metallicity(field, data):
                ret = (
                    data[(ptype, "ElevenMetalMasses")][:, 1].sum(axis=1)
                    / data[(ptype, "particle_mass")]
                )
                return ret

            def _h_fraction(field, data):
                ret = (
                    data[(ptype, "ElevenMetalMasses")].sum(axis=1)
                    / data[(ptype, "particle_mass")]
                )
                return 1.0 - ret

            elem_names = elem_names_opts[11]

        elif fname == "MetalMasses":
            n_elem = len(
                [
                    fd
                    for fd in self.ds.field_list
                    if fd[0] == ptype and fd[1].startswith("MetalMasses")
                ]
            )
            elem_names = elem_names_opts[n_elem]
            no_He = "He" not in elem_names

            def _fraction(field, data, i: int):
                return (
                    data[(ptype, f"MetalMasses_{i:02d}")]
                    / data[(ptype, "particle_mass")]
                )

            def _metallicity(field, data):
                mass = 0.0
                start_idx = int(not no_He)
                for i in range(start_idx, n_elem):
                    mass += data[(ptype, f"MetalMasses_{i:02d}")]
                ret = mass / data[(ptype, "particle_mass")]
                return ret

            if no_He:
                _h_fraction = None

            else:

                def _h_fraction(field, data):
                    mass = 0.0
                    for i in range(n_elem):
                        mass += data[(ptype, f"MetalMasses_{i:02d}")]
                    ret = mass / data[(ptype, "particle_mass")]
                    return 1.0 - ret

        else:
            raise KeyError(
                f"Making element fraction fields from '{ptype}','{fname}' not possible!"
            )
        return _fraction, _h_fraction, _metallicity, elem_names

    def _setup_metal_masses(self, ptype, fname):
        """
        This function breaks the ElevenMetalMasses field (if present)
        into its eleven component metal fraction fields and adds
        corresponding metal density fields which will later get smoothed

        This gets used with the magneticum_box2_hr format
        as defined in the gadget_field_specs in frontends/gadget/definitions.py
        """
        sampling_type = "local" if ptype in self.ds._sph_ptypes else "particle"
        (
            _fraction,
            _h_fraction,
            _metallicity,
            elem_names,
        ) = self._make_fraction_functions(ptype, fname)

        def _metal_density(field, data, elem_name: str):
            return data[(ptype, f"{elem_name}_fraction")] * data[(ptype, "density")]

        for i, elem_name in enumerate(elem_names):
            # add the element fraction fields
            self.add_field(
                (ptype, f"{elem_name}_fraction"),
                sampling_type=sampling_type,
                function=partial(_fraction, i=i),
                units="",
            )

            # add the element density fields
            self.add_field(
                (ptype, f"{elem_name}_density"),
                sampling_type=sampling_type,
                function=partial(_metal_density, elem_name=elem_name),
                units=self.ds.unit_system["density"],
            )

        # metallicity
        self.add_field(
            (ptype, "metallicity"),
            sampling_type=sampling_type,
            function=_metallicity,
            units="",
        )

        if _h_fraction is None:
            # no helium, so can't compute hydrogen
            species_names = elem_names[-1]

        else:
            # hydrogen fraction and density
            self.add_field(
                (ptype, "H_fraction"),
                sampling_type=sampling_type,
                function=_h_fraction,
                units="",
            )

            def _h_density(field, data):
                return data[(ptype, "H_fraction")] * data[(ptype, "density")]

            self.add_field(
                (ptype, "H_density"),
                sampling_type=sampling_type,
                function=_h_density,
                units=self.ds.unit_system["density"],
            )

            species_names = ["H"] + elem_names[:-1]

        if "Ej" in elem_names:

            def _ej_mass(field, data):
                return data[(ptype, "Ej_fraction")] * data[(ptype, "particle_mass")]

            self.add_field(
                (ptype, "Ej_mass"),
                sampling_type=sampling_type,
                function=_ej_mass,
                units=self.ds.unit_system["mass"],
            )
            if sampling_type == "local":
                self.alias(("gas", "Ej_mass"), (ptype, "Ej_mass"))

        return species_names

    def _check_whitelist_species_fields(self, ptype):
        species_names = []
        for field in self.ds.field_list:
            if (
                field[0] == ptype
                and field[1].endswith(("_fraction", "_density"))
                and field[1] in sph_whitelist_fields
            ):
                symbol, _, _ = field[1].partition("_")
                species_names.append(symbol)
        return sorted(species_names, key=lambda symbol: periodic_table[symbol].num)

    def setup_gas_particle_fields(self, ptype):
        if (ptype, "Temperature") not in self.ds.field_list:
            if (ptype, "ElectronAbundance") in self.ds.field_list:

                def _temperature(field, data):
                    # Assume cosmic abundances
                    x_H = _primordial_mass_fraction["H"]
                    gamma = 5.0 / 3.0
                    a_e = data[ptype, "ElectronAbundance"]
                    mu = 4.0 / (3.0 * x_H + 1.0 + 4.0 * x_H * a_e)
                    ret = data[ptype, "InternalEnergy"] * (gamma - 1) * mu * mp / kb
                    return ret.in_units(self.ds.unit_system["temperature"])

            else:

                def _temperature(field, data):
                    gamma = 5.0 / 3.0
                    ret = (
                        data[ptype, "InternalEnergy"]
                        * (gamma - 1)
                        * data.ds.mu
                        * mp
                        / kb
                    )
                    return ret.in_units(self.ds.unit_system["temperature"])

            self.add_field(
                (ptype, "Temperature"),
                sampling_type="particle",
                function=_temperature,
                units=self.ds.unit_system["temperature"],
            )

        self.alias((ptype, "temperature"), (ptype, "Temperature"))
        # need to do this manually since that automatic aliasing that happens
        # in the FieldInfoContainer base class has already happened at this
        # point
        self.alias(("gas", "temperature"), (ptype, "Temperature"))
