"""
Gizmo-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.magnetic_field import \
    setup_magnetic_field_aliases
from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field
from yt.fields.species_fields import \
    add_species_field_by_density, \
    setup_species_fields
from yt.frontends.gadget.fields import \
    GadgetFieldInfo
from yt.frontends.sph.fields import \
    SPHFieldInfo

metal_elements = ["He", "C", "N", "O", "Ne",
                  "Mg", "Si", "S", "Ca", "Fe"]

class GizmoFieldInfo(GadgetFieldInfo):
    # The known fields list is according to the GIZMO User Guide. See
    # http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html#snaps-reading
    known_particle_fields = (
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("InternalEnergy", ("code_specific_energy", ["thermal_energy"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("ElectronAbundance", ("", [], None)),
        ("NeutralHydrogenAbundance", ("", [], None)),
        ("StarFormationRate", ("Msun / yr", [], None)),
        ("Metallicity", ("code_metallicity", ["metallicity"], None)),
        ("Metallicity_00", ("", ["metallicity"], None)),
        ("Metallicity_01", ("", ["He_metallicity"], None)),
        ("Metallicity_02", ("", ["C_metallicity"], None)),
        ("Metallicity_03", ("", ["N_metallicity"], None)),
        ("Metallicity_04", ("", ["O_metallicity"], None)),
        ("Metallicity_05", ("", ["Ne_metallicity"], None)),
        ("Metallicity_06", ("", ["Mg_metallicity"], None)),
        ("Metallicity_07", ("", ["Si_metallicity"], None)),
        ("Metallicity_08", ("", ["S_metallicity"], None)),
        ("Metallicity_09", ("", ["Ca_metallicity"], None)),
        ("Metallicity_10", ("", ["Fe_metallicity"], None)),
        ("ArtificialViscosity", ("", [], None)),
        ("MagneticField",
         ("code_magnetic", ["particle_magnetic_field"], None)),
        ("DivergenceOfMagneticField",
         ("code_magnetic / code_length", [], None)),
        ("StellarFormationTime", ("", [], None)),
        # "StellarFormationTime" has different meanings in (non-)cosmological
        # runs, so units are left blank here.
        ("BH_Mass", ("code_mass", [], None)),
        ("BH_Mdot", ("code_mass / code_time", [], None)),
        ("BH_Mass_AlphaDisk", ("code_mass", [], None)),
    )

    def __init__(self, *args, **kwargs):
        super(SPHFieldInfo, self).__init__(*args, **kwargs)
        if ("PartType0", "Metallicity_00") in self.field_list:
            self.nuclei_names = metal_elements
            self.species_names = ["H", "H_p1"] + metal_elements

    def setup_particle_fields(self, ptype):
        FieldInfoContainer.setup_particle_fields(self, ptype)
        if ptype in ("PartType0",):
            self.setup_gas_particle_fields(ptype)
            setup_species_fields(self, ptype)
        if ptype in ("PartType4",):
            self.setup_star_particle_fields(ptype)

    def setup_gas_particle_fields(self, ptype):
        super(GizmoFieldInfo, self).setup_gas_particle_fields(ptype)

        def _h_density(field, data):
            x_H = 1.0 - data[(ptype, "He_metallicity")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              data[(ptype, "NeutralHydrogenAbundance")]

        self.add_field(
            (ptype, "H_density"),
            sampling_type="particle",
            function=_h_density,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H", particle_type=True)
        for suffix in ["density", "fraction", "mass", "number_density"]:
            self.alias((ptype, "H_p0_%s" % suffix), (ptype, "H_%s" % suffix))

        def _h_p1_density(field, data):
            x_H = 1.0 - data[(ptype, "He_metallicity")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              (1.0 - data[(ptype, "NeutralHydrogenAbundance")])

        self.add_field(
            (ptype, "H_p1_density"),
            sampling_type="particle",
            function=_h_p1_density,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H_p1", particle_type=True)

        def _nuclei_mass_density_field(field, data):
            species = field.name[1][:field.name[1].find("_")]
            return data[ptype, "density"] * \
              data[ptype, "%s_metallicity" % species]

        num_neighbors = 64
        for species in ['H', 'H_p0', 'H_p1']:
            for suf in ["_density", "_number_density"]:
                field = "%s%s" % (species, suf)
                fn = add_volume_weighted_smoothed_field(
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", field,
                    self, num_neighbors)
                self.alias(("gas", field), fn[0])

        for species in self.nuclei_names:
            self.add_field(
                (ptype, "%s_nuclei_mass_density" % species),
                sampling_type="particle",
                function=_nuclei_mass_density_field,
                units=self.ds.unit_system["density"])

            for suf in ["_nuclei_mass_density", "_metallicity"]:
                field = "%s%s" % (species, suf)
                fn = add_volume_weighted_smoothed_field(
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", field,
                    self, num_neighbors)

                self.alias(("gas", field), fn[0])

        magnetic_field = "MagneticField"
        if (ptype, magnetic_field) in self.field_list:
            setup_magnetic_field_aliases(
                self, ptype, magnetic_field, ftype=ptype
            )

    def setup_star_particle_fields(self, ptype):
        def _creation_time(field, data):
            if data.ds.cosmological_simulation:
                a_form = data['StellarFormationTime']
                z_form = 1 / a_form - 1
                creation_time = data.ds.cosmology.t_from_z(z_form)
            else:
                t_form = data["StellarFormationTime"]
                creation_time = data.ds.arr(t_form, "code_time")
            return creation_time

        self.add_field(
            (ptype, "creation_time"),
            sampling_type="particle",
            function=_creation_time,
            units=self.ds.unit_system["time"])

        def _age(field, data):
            return data.ds.current_time - data["creation_time"]

        self.add_field(
            (ptype, "age"),
            sampling_type="particle",
            function=_age,
            units=self.ds.unit_system["time"])
