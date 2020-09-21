"""
Skeleton-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
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

class CM1FieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("uinterp", ("m/s",    ["velocity_x"], None)),
        ("vinterp", ("m/s",    ["velocity_y"], None)),
        ("winterp", ("m/s",    ["velocity_z"], None)),
        ("hwin_sr", ("m/s",    ["storm_relative_horizontal_wind_speed"], None)),
        ("windmag_sr", ("m/s", ["storm_relative_3D_wind_speed"], None)),
        ("hwin_gr", ("m/s",    ["ground_relative_horizontal_wind_speed"], None)),
        ("thpert", ("K",       ["potential_temperature_perturbation"], None)),
        ("thrhopert", ("K",    ["density_potential_temperature_perturbation"], None)),
        ("prespert", ("hPa",   ["presure_perturbation"], None)),
        ("rhopert", ("kg/m**3", ["density_perturbation"], None)),
        ("dbz", ("dBZ",        ["simulated_reflectivity"], None)),
        ("qvpert", ("g/kg",    ["water_vapor_mixing_ratio_perturbation"], None)),
        ("qc", ("g/kg",        ["cloud_liquid_water_mixing_ratio"], None)),
        ("qr", ("g/kg",        ["rain_mixing_ratio"], None)),
        ("qi", ("g/kg",        ["cloud_ice_mixing_ratio"], None)),
        ("qs", ("g/kg",        ["snow_mixing_ratio"], None)),
        ("qg", ("g/kg",        ["graupel_or_hail_mixing_ratio"], None)),
        ("qcloud", ("g/kg",    ["sum_of_cloud_water_and_cloud_ice_mixing_ratios"], None)),
        ("qprecip", ("g/kg",   ["sum_of_rain_graupel_snow_mixing_ratios"], None)),
        ("nci", ("1/cm**3",      ["number_concerntration_of_cloud_ice"], None)),
        ("ncr", ("1/cm**3",      ["number_concentration_of_rain"], None)),
        ("ncs", ("1/cm**3",      ["number_concentration_of_snow"], None)),
        ("ncg", ("1/cm**3",      ["number_concentration_of_graupel_or_hail"], None)),
        ("xvort", ("1/s",     ["vorticity_x"], None)),
        ("yvort", ("1/s",     ["vorticity_y"], None)),
        ("zvort", ("1/s",     ["vorticity_z"], None)),
        ("hvort", ("1/s",     ["horizontal_vorticity_magnitude"], None)),
        ("vortmag", ("1/s",   ["vorticity_magnitude"], None)),
        ("streamvort", ("1/s",["streamwise_vorticity"], None)),
        ("khh", ("m**2/s",      ["khh"], None)),
        ("khv", ("m**2/s",      ["khv"], None)),
        ("kmh", ("m**2/s",      ["kmh"], None)),
        ("kmv", ("m**2/s",      ["kmv"], None))
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def __init__(self, ds, field_list):
        super(CM1FieldInfo, self).__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        pass

    def setup_particle_fields(self, ptype):
        super(CM1FieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
