from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.fields.species_fields import add_species_field_by_fraction, setup_species_fields
from yt.frontends.gadget.api import GadgetFieldInfo
from yt.utilities.chemical_formulas import ChemicalFormula
from yt.utilities.physical_ratios import _primordial_mass_fraction

metal_elements = ["He", "C", "N", "O", "Ne", "Mg", "Si", "Fe"]


class ArepoFieldInfo(GadgetFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        if ds.cosmological_simulation:
            GFM_SFT_units = "dimensionless"
        else:
            GFM_SFT_units = "code_length/code_velocity"
        self.known_particle_fields += (
            ("GFM_StellarFormationTime", (GFM_SFT_units, ["stellar_age"], None)),
            ("MagneticField", ("code_magnetic", ["particle_magnetic_field"], None)),
            (
                "MagneticFieldDivergence",
                ("code_magnetic/code_length", ["magnetic_field_divergence"], None),
            ),
            ("GFM_Metallicity", ("", ["metallicity"], None)),
            ("GFM_Metals_00", ("", ["H_fraction"], None)),
            ("GFM_Metals_01", ("", ["He_fraction"], None)),
            ("GFM_Metals_02", ("", ["C_fraction"], None)),
            ("GFM_Metals_03", ("", ["N_fraction"], None)),
            ("GFM_Metals_04", ("", ["O_fraction"], None)),
            ("GFM_Metals_05", ("", ["Ne_fraction"], None)),
            ("GFM_Metals_06", ("", ["Mg_fraction"], None)),
            ("GFM_Metals_07", ("", ["Si_fraction"], None)),
            ("GFM_Metals_08", ("", ["Fe_fraction"], None)),
            ("GFM_StellarPhotometrics_00", ("", ["U_magnitude"], None)),
            ("GFM_StellarPhotometrics_01", ("", ["B_magnitude"], None)),
            ("GFM_StellarPhotometrics_02", ("", ["V_magnitude"], None)),
            ("GFM_StellarPhotometrics_03", ("", ["K_magnitude"], None)),
            ("GFM_StellarPhotometrics_04", ("", ["g_magnitude"], None)),
            ("GFM_StellarPhotometrics_05", ("", ["r_magnitude"], None)),
            ("GFM_StellarPhotometrics_06", ("", ["i_magnitude"], None)),
            ("GFM_StellarPhotometrics_07", ("", ["z_magnitude"], None)),
            (
                "CosmicRaySpecificEnergy",
                ("code_specific_energy", ["specific_cosmic_ray_energy"], None),
            ),
        )
        super().__init__(ds, field_list, slice_info=slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):
        FieldInfoContainer.setup_particle_fields(self, ptype)
        if ptype == "PartType0":
            self.setup_gas_particle_fields(ptype)
            setup_species_fields(self, ptype)

    def setup_gas_particle_fields(self, ptype):
        super().setup_gas_particle_fields(ptype)

        # Since the AREPO gas "particles" are Voronoi cells, we can
        # define a volume here
        def _volume(field, data):
            return data[ptype, "mass"] / data[ptype, "density"]

        self.add_field(
            (ptype, "cell_volume"),
            function=_volume,
            sampling_type="local",
            units=self.ds.unit_system["volume"],
        )

        if (ptype, "InternalEnergy") in self.field_list:

            def _pressure(field, data):
                return (
                    (data.ds.gamma - 1.0)
                    * data[ptype, "density"]
                    * data[ptype, "InternalEnergy"]
                )

            self.add_field(
                (ptype, "pressure"),
                function=_pressure,
                sampling_type="particle",
                units=self.ds.unit_system["pressure"],
            )

            self.alias((ptype, "pressure"), ("gas", "pressure"))

        if (ptype, "GFM_Metals_00") in self.field_list:
            self.nuclei_names = metal_elements
            self.species_names = ["H"] + metal_elements

        if (ptype, "MagneticField") in self.field_list:
            setup_magnetic_field_aliases(self, ptype, "MagneticField")

        if (ptype, "NeutralHydrogenAbundance") in self.field_list:

            def _h_p0_fraction(field, data):
                return (
                    data[ptype, "GFM_Metals_00"]
                    * data[ptype, "NeutralHydrogenAbundance"]
                )

            self.add_field(
                (ptype, "H_p0_fraction"),
                sampling_type="particle",
                function=_h_p0_fraction,
                units="",
            )

            def _h_p1_fraction(field, data):
                return data[ptype, "GFM_Metals_00"] * (
                    1.0 - data[ptype, "NeutralHydrogenAbundance"]
                )

            self.add_field(
                (ptype, "H_p1_fraction"),
                sampling_type="particle",
                function=_h_p1_fraction,
                units="",
            )

            add_species_field_by_fraction(self, ptype, "H_p0")
            add_species_field_by_fraction(self, ptype, "H_p1")

            for species in ["H", "H_p0", "H_p1"]:
                for suf in ["_density", "_number_density"]:
                    field = f"{species}{suf}"
                    self.alias(("gas", field), (ptype, field))

        if (ptype, "ElectronAbundance") in self.field_list:

            # If we have ElectronAbundance but not NeutralHydrogenAbundance, assume the
            # cosmic value for hydrogen to generate the H_number_density
            if (ptype, "NeutralHydrogenAbundance") not in self.field_list:
                amu_cgs = self.ds.units.physical_constants.amu_cgs
                muinv = _primordial_mass_fraction["H"] / ChemicalFormula("H").weight

                def _h_number_density(field, data):
                    return data["gas", "density"] * muinv / amu_cgs

                self.add_field(
                    (ptype, "H_number_density"),
                    sampling_type="particle",
                    function=_h_number_density,
                    units=self.ds.unit_system["number_density"],
                )
                self.alias(("gas", "H_number_density"), (ptype, "H_number_density"))
                self.alias(("gas", "H_nuclei_density"), ("gas", "H_number_density"))

            def _el_number_density(field, data):
                return (
                    data[ptype, "ElectronAbundance"] * data[ptype, "H_number_density"]
                )

            self.add_field(
                (ptype, "El_number_density"),
                sampling_type="particle",
                function=_el_number_density,
                units=self.ds.unit_system["number_density"],
            )
            self.alias(("gas", "El_number_density"), (ptype, "El_number_density"))

        if (ptype, "CosmicRaySpecificEnergy") in self.field_list:

            self.alias(
                (ptype, "specific_cosmic_ray_energy"),
                ("gas", "specific_cosmic_ray_energy"),
            )

            def _cr_energy_density(field, data):
                return (
                    data["PartType0", "specific_cosmic_ray_energy"]
                    * data["gas", "density"]
                )

            self.add_field(
                ("gas", "cosmic_ray_energy_density"),
                _cr_energy_density,
                sampling_type="local",
                units=self.ds.unit_system["pressure"],
            )

            def _cr_pressure(field, data):
                return (data.ds.gamma_cr - 1.0) * data[
                    "gas", "cosmic_ray_energy_density"
                ]

            self.add_field(
                ("gas", "cosmic_ray_pressure"),
                _cr_pressure,
                sampling_type="local",
                units=self.ds.unit_system["pressure"],
            )
