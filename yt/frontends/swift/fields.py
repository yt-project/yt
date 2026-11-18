from functools import partial

from yt.frontends.sph.fields import SPHFieldInfo

# Column index of each species inside the packed SWIFT "SpeciesFractions"
# particle field, as produced by the SWIFT chemistry/cooling network.
# Elements map to a list of (ion_name, column_index); molecules map
# straight to a column index.
ELEMENT_MAP = {
    "H": [("HI", 1), ("HII", 2), ("Hm", 3)],
    "He": [("HeI", 4), ("HeII", 5), ("HeIII", 6)],
    "C": [
        ("CI", 7),
        ("CII", 8),
        ("CIII", 9),
        ("CIV", 10),
        ("CV", 11),
        ("CVI", 12),
        ("CVII", 13),
        ("Cm", 14),
    ],
    "N": [
        ("NI", 15),
        ("NII", 16),
        ("NIII", 17),
        ("NIV", 18),
        ("NV", 19),
        ("NVI", 20),
        ("NVII", 21),
        ("NVIII", 22),
    ],
    "O": [
        ("OI", 23),
        ("OII", 24),
        ("OIII", 25),
        ("OIV", 26),
        ("OV", 27),
        ("OVI", 28),
        ("OVII", 29),
        ("OVIII", 30),
        ("OIX", 31),
        ("Om", 32),
    ],
    "Ne": [
        ("NeI", 33),
        ("NeII", 34),
        ("NeIII", 35),
        ("NeIV", 36),
        ("NeV", 37),
        ("NeVI", 38),
        ("NeVII", 39),
        ("NeVIII", 40),
        ("NeIX", 41),
        ("NeX", 42),
        ("NeXI", 43),
    ],
    "Mg": [
        ("MgI", 44),
        ("MgII", 45),
        ("MgIII", 46),
        ("MgIV", 47),
        ("MgV", 48),
        ("MgVI", 49),
        ("MgVII", 50),
        ("MgVIII", 51),
        ("MgIX", 52),
        ("MgX", 53),
        ("MgXI", 54),
        ("MgXII", 55),
        ("MgXIII", 56),
    ],
    "Si": [
        ("SiI", 57),
        ("SiII", 58),
        ("SiIII", 59),
        ("SiIV", 60),
        ("SiV", 61),
        ("SiVI", 62),
        ("SiVII", 63),
        ("SiVIII", 64),
        ("SiIX", 65),
        ("SiX", 66),
        ("SiXI", 67),
        ("SiXII", 68),
        ("SiXIII", 69),
        ("SiXIV", 70),
        ("SiXV", 71),
    ],
    "S": [
        ("SI", 72),
        ("SII", 73),
        ("SIII", 74),
        ("SIV", 75),
        ("SV", 76),
        ("SVI", 77),
        ("SVII", 78),
        ("SVIII", 79),
        ("SIX", 80),
        ("SX", 81),
        ("SXI", 82),
        ("SXII", 83),
        ("SXIII", 84),
        ("SXIV", 85),
        ("SXV", 86),
        ("SXVI", 87),
        ("SXVII", 88),
    ],
    "Ca": [
        ("CaI", 89),
        ("CaII", 90),
        ("CaIII", 91),
        ("CaIV", 92),
        ("CaV", 93),
        ("CaVI", 94),
        ("CaVII", 95),
        ("CaVIII", 96),
        ("CaIX", 97),
        ("CaX", 98),
        ("CaXI", 99),
        ("CaXII", 100),
        ("CaXIII", 101),
        ("CaXIV", 102),
        ("CaXV", 103),
        ("CaXVI", 104),
        ("CaXVII", 105),
        ("CaXVIII", 106),
        ("CaXIX", 107),
        ("CaXX", 108),
        ("CaXXI", 109),
    ],
    "Fe": [
        ("FeI", 110),
        ("FeII", 111),
        ("FeIII", 112),
        ("FeIV", 113),
        ("FeV", 114),
        ("FeVI", 115),
        ("FeVII", 116),
        ("FeVIII", 117),
        ("FeIX", 118),
        ("FeX", 119),
        ("FeXI", 120),
        ("FeXII", 121),
        ("FeXIII", 122),
        ("FeXIV", 123),
        ("FeXV", 124),
        ("FeXVI", 125),
        ("FeXVII", 126),
        ("FeXVIII", 127),
        ("FeXIX", 128),
        ("FeXX", 129),
        ("FeXXI", 130),
        ("FeXXII", 131),
        ("FeXXIII", 132),
        ("FeXXIV", 133),
        ("FeXXV", 134),
        ("FeXXVI", 135),
        ("FeXXVII", 136),
    ],
}

MOLECULE_MAP = {
    "H2": 137,
    "H2p": 138,
    "H3p": 139,
    "OH": 140,
    "H2O": 141,
    "C2": 142,
    "O2": 143,
    "HCOp": 144,
    "CH": 145,
    "CH2": 146,
    "CH3p": 147,
    "CO": 148,
    "CHp": 149,
    "CH2p": 150,
    "OHp": 151,
    "H2Op": 152,
    "H3Op": 153,
    "COp": 154,
    "HOCp": 155,
    "O2p": 156,
}


class SwiftFieldInfo(SPHFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        self.known_particle_fields += (
            (
                "InternalEnergies",
                ("code_specific_energy", ["specific_thermal_energy"], None),
            ),
            ("Densities", ("code_mass / code_length**3", ["density"], None)),
            ("SmoothingLengths", ("code_length", ["smoothing_length"], None)),
        )
        super().__init__(ds, field_list, slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super().setup_particle_fields(ptype, *args, **kwargs)

        if ptype in ("PartType0", "Gas"):
            self.setup_gas_particle_fields(ptype)
            if (ptype, "SpeciesFractions") in self.ds.field_list:
                self.species_names = self._setup_species_fractions(ptype)

    def _setup_species_fractions(self, ptype):
        """
        Splits the packed "SpeciesFractions" particle field (one column per
        ionization state / molecule, as produced by SWIFT's chemistry
        network) into two fields per ion/molecule:

        - "<species>_fraction_raw": the unmodified column value.
        - "<species>_fraction": for ions, this is renormalized so it sums
          to 1 across all ionization states of the parent element (i.e.
          "fraction of carbon atoms that are C IV"). For molecules, which
          have no parent-element group to normalize against, this is
          identical to "<species>_fraction_raw".
        """

        def _raw_value(field, data, idx: int):
            return data[ptype, "SpeciesFractions"][:, idx]

        def _ion_fraction(field, data, idx: int, indices: list):
            all_species = data[ptype, "SpeciesFractions"]
            num = all_species[:, idx]
            total = all_species[:, indices].sum(axis=1)
            out = num.copy()
            mask = total > 0
            out[mask] = num[mask] / total[mask]
            out[~mask] = 0
            return out

        species_names = []

        for _element, ions in ELEMENT_MAP.items():
            indices = [idx for _, idx in ions]
            for name, idx in ions:
                self.add_field(
                    (ptype, f"{name}_fraction_raw"),
                    sampling_type="particle",
                    function=partial(_raw_value, idx=idx),
                    units="",
                )
                self.add_field(
                    (ptype, f"{name}_fraction"),
                    sampling_type="particle",
                    function=partial(_ion_fraction, idx=idx, indices=indices),
                    units="",
                )
                species_names.append(name)

        for name, idx in MOLECULE_MAP.items():
            self.add_field(
                (ptype, f"{name}_fraction_raw"),
                sampling_type="particle",
                function=partial(_raw_value, idx=idx),
                units="",
            )
            self.alias((ptype, f"{name}_fraction"), (ptype, f"{name}_fraction_raw"))
            species_names.append(name)

        return species_names

    def setup_gas_particle_fields(self, ptype):
        self.alias((ptype, "temperature"), (ptype, "Temperatures"))
        self.alias(("gas", "temperature"), (ptype, "Temperatures"))

        for ax in ("x", "y", "z"):
            self.alias((ptype, ax), (ptype, "particle_position_" + ax))
