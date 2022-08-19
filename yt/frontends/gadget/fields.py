from functools import partial

from yt.fields.particle_fields import sph_whitelist_fields
from yt.frontends.sph.fields import SPHFieldInfo
from yt.utilities.periodic_table import periodic_table
from yt.utilities.physical_constants import kb, mp
from yt.utilities.physical_ratios import _primordial_mass_fraction


class GadgetFieldInfo(SPHFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        if ds.gen_hsmls:
            hsml = (("smoothing_length", ("code_length", [], None)),)
            self.known_particle_fields += hsml
        super().__init__(ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):

        # setup some special fields that only make sense for SPH particles

        if (ptype, "FourMetalFractions") in self.ds.field_list:
            self.species_names = self._setup_four_metal_fractions(ptype)
        elif (ptype, "ElevenMetalMasses") in self.ds.field_list:
            self.species_names = self._setup_eleven_metal_masses(ptype)
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
        metal_names = ["C", "O", "Si", "Fe"]

        def _Fraction(field, data, i: int):
            return data[(ptype, "FourMetalFractions")][:, i]

        def _Metal_density(field, data, i: int):
            return data[(ptype, "FourMetalFractions")][:, i] * data[(ptype, "density")]

        for i, metal_name in enumerate(metal_names):

            # add the metal fraction fields
            self.add_field(
                (ptype, metal_name + "_fraction"),
                sampling_type="particle",
                function=partial(_Fraction, i=i),
                units="",
            )

            # add the metal density fields
            self.add_field(
                (ptype, metal_name + "_density"),
                sampling_type="particle",
                function=partial(_Metal_density, i=i),
                units=self.ds.unit_system["density"],
            )

        return metal_names

    def _setup_eleven_metal_masses(self, ptype):
        """
        This function breaks the ElevenMetalMasses field (if present)
        into its eleven component metal fraction fields and adds
        corresponding metal density fields which will later get smoothed

        This gets used with the magneticum_box2_hr format
        as defined in the gadget_field_specs in frontends/gadget/definitions.py
        """
        metal_names = ["He", "C", "Ca", "O", "N", "Ne", "Mg", "S", "Si", "Fe", "Ej"]

        def _Fraction(field, data, i: int):
            return data[(ptype, "ElevenMetalMasses")][:, i] / data[(ptype, "Mass")]

        def _Metal_density(field, data, metal_name: str):
            return data[(ptype, metal_name + "_fraction")] * data[(ptype, "density")]

        for i, metal_name in enumerate(metal_names):

            # add the metal fraction fields
            self.add_field(
                (ptype, metal_name + "_fraction"),
                sampling_type="particle",
                function=partial(_Fraction, i=i),
                units="",
            )

            # add the metal density fields
            self.add_field(
                (ptype, metal_name + "_density"),
                sampling_type="particle",
                function=partial(_Metal_density, metal_name=metal_name),
                units=self.ds.unit_system["density"],
            )

        # hydrogen fraction and density
        def _h_fraction(field, data):
            ret = data[(ptype, "ElevenMetalMasses")].sum(axis=1) / data[(ptype, "Mass")]
            return 1.0 - ret

        self.add_field(
            (ptype, "H_fraction"),
            sampling_type="particle",
            function=_h_fraction,
            units="",
        )

        def _h_density(field, data):
            return data[(ptype, "H_fraction")] * data[(ptype, "density")]

        self.add_field(
            (ptype, "H_density"),
            sampling_type="particle",
            function=_h_density,
            units=self.ds.unit_system["density"],
        )

        return ["H"] + metal_names[:-1]

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
