from yt.frontends.sph.fields import SPHFieldInfo
from yt.utilities.physical_constants import kb, mp
from yt.utilities.physical_ratios import _primordial_mass_fraction


class GadgetFieldInfo(SPHFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        if ds.gen_hsmls:
            hsml = (("smoothing_length", ("code_length", [], None)),)
            self.known_particle_fields += hsml
        super(GadgetFieldInfo, self).__init__(ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):

        # setup some special fields that only make sense for SPH particles

        if (ptype, "FourMetalFractions") in self.ds.field_list:
            self._setup_four_metal_fractions(ptype)

        super(GadgetFieldInfo, self).setup_particle_fields(ptype, *args, **kwargs)

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
        for i, metal_name in enumerate(metal_names):

            # add the metal fraction fields
            def _Fraction_wrap(i):
                def _Fraction(field, data):
                    return data[(ptype, "FourMetalFractions")][:, i]

                return _Fraction

            self.add_field(
                (ptype, metal_name + "_fraction"),
                sampling_type="particle",
                function=_Fraction_wrap(i),
                units="",
            )

            # add the metal density fields
            def _Density_wrap(i):
                def _Metal_density(field, data):
                    return (
                        data[(ptype, "FourMetalFractions")][:, i]
                        * data[(ptype, "density")]
                    )

                return _Metal_density

            self.add_field(
                (ptype, metal_name + "_density"),
                sampling_type="particle",
                function=_Density_wrap(i),
                units=self.ds.unit_system["density"],
            )

    def setup_gas_particle_fields(self, ptype):
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
                ret = data[ptype, "InternalEnergy"] * (gamma - 1) * data.ds.mu * mp / kb
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
