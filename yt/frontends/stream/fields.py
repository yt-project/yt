"""
Fields specific to Streaming data



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
import yt.fields.api
from yt.fields.particle_fields import \
    particle_deposition_functions, \
    particle_vector_functions

class StreamFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("code_mass/code_length**3", [], None)),
        ("dark_matter_density", ("code_mass/code_length**3", [], None)),
        ("number_density", ("1/code_length**3", [], None)),
        ("pressure", ("dyne/code_length**2", [], None)),
        ("thermal_energy", ("erg / g", [], None)),
        ("temperature", ("K", [], None)),
        ("velocity_x", ("code_length/code_time", [], None)),
        ("velocity_y", ("code_length/code_time", [], None)),
        ("velocity_z", ("code_length/code_time", [], None)),
        ("magnetic_field_x", ("gauss", [], None)),
        ("magnetic_field_y", ("gauss", [], None)),
        ("magnetic_field_z", ("gauss", [], None)),
        ("radiation_acceleration_x", ("code_length/code_time**2", [], None)),
        ("radiation_acceleration_y", ("code_length/code_time**2", [], None)),
        ("radiation_acceleration_z", ("code_length/code_time**2", [], None)),

        # We need to have a bunch of species fields here, too
        ("metal_density",   ("code_mass/code_length**3", [], None)),
        ("hi_density",      ("code_mass/code_length**3", [], None)),
        ("hii_density",     ("code_mass/code_length**3", [], None)),
        ("h2i_density",     ("code_mass/code_length**3", [], None)),
        ("h2ii_density",    ("code_mass/code_length**3", [], None)),
        ("h2m_density",     ("code_mass/code_length**3", [], None)),
        ("hei_density",     ("code_mass/code_length**3", [], None)),
        ("heii_density",    ("code_mass/code_length**3", [], None)),
        ("heiii_density",   ("code_mass/code_length**3", [], None)),
        ("hdi_density",     ("code_mass/code_length**3", [], None)),
        ("di_density",      ("code_mass/code_length**3", [], None)),
        ("dii_density",     ("code_mass/code_length**3", [], None)),
    )

    known_particle_fields = (
        ("particle_position", ("code_length", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity", ("code_length/code_time", [], None)),
        ("particle_velocity_x", ("code_length/code_time", [], None)),
        ("particle_velocity_y", ("code_length/code_time", [], None)),
        ("particle_velocity_z", ("code_length/code_time", [], None)),
        ("particle_index", ("", [], None)),
        ("particle_gas_density", ("code_mass/code_length**3", [], None)),
        ("particle_gas_temperature", ("K", [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("smoothing_length", ("code_length", [], None)),
        ("density", ("code_mass/code_length**3", [], None)),
        ("temperature", ("code_temperature", [], None)),
        ("creation_time", ("code_time", [], None)),
    )

    def setup_fluid_fields(self):
        for field in self.ds.stream_handler.field_units:
            units = self.ds.stream_handler.field_units[field]
            if units != '': self.add_output_field(field, units=units)

    def add_output_field(self, name, **kwargs):
        if name in self.ds.stream_handler.field_units:
            kwargs['units'] = self.ds.stream_handler.field_units[name]
        super(StreamFieldInfo, self).add_output_field(name, **kwargs)
