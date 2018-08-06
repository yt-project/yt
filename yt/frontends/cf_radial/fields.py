"""
CF-radial-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class CFRadialFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("reflectivity", ("dBZ", [], None)),
        ("correct_reflectivity", ("dBZ", [], None)),
        ("total_power", ("dBZ", [], None)),
        ("velocity", (
            "m/s",
            ["mean_doppler_velocity", "doppler_velocity"], None)),
        ("corrected_velocity", ("m/s", [], None)),
        ("simulated_velocity", ("m/s", [], None)),
        ("spectrum_width", ("m/s", [], None)),
        ("differential_reflectivity", ("dB", [], None)),
        ("corrected_differential_reflectivity", ("dB", [], None)),
        ("cross_correlation_ratio", ("ratio", [], None)),
        ("normalized_coherent_power", ("ratio", [], None)),
        ("differential_phase", ("degree", [], None)),
        ("unfolded_differential_phase", ("degree", [], None)),
        ("corrected_differential_phase", ("degree", [], None)),
        ("specific_differential_phase", ("degree/km", [], None)),
        ("corrected_specific_differential_phase", ("degree/km", [], None)),
        ("linear_depolarization_ratio", ("dB", [], None)),
        ("linear_depolarization_ratio_h", ("dB", [], None)),
        ("linear_depolarization_ratio_v", ("dB", [], None)),
        ("signal_to_noise_ratio", ("dB", [], None)),
        ("rain_rate", ("kg/m2/s", [], None)),
        ("radar_estimated_rain_rate", ("mm/hr", [], None)),
        ("radar_echo_classification", ("legend", [], None)),
        ("specific_attenuation", ("dB/km", [], None)),
        ("differential_phase_texture", ("degree", [], None)),
        ("eastward_wind_component", ("m/s", [], None)),
        ("northward_wind_component", ("m/s", [], None)),
        ("vertical_wind_component", ("m/s", [], None))
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def __init__(self, ds, field_list):
        super(CFRadialFieldInfo, self).__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        pass

    def setup_particle_fields(self, ptype):
        super(CFRadialFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
