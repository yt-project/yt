"""
OWLS-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import *
from yt.fields.field_info_container import \
    FieldInfoContainer
from .definitions import \
    gadget_ptypes, \
    ghdf5_ptypes

from yt.fields.species_fields import add_species_field_by_fraction



# Here are helper functions for things like vector fields and so on.

def _get_conv(cf):
    def _convert(data):
        return data.convert(cf)
    return _convert

class SPHFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ("Mass", ("code_mass", ["particle_mass"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocity", ("code_velocity", ["particle_velocity"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("InternalEnergy", ("", ["thermal_energy"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("MaximumTemperature", ("K", [], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Epsilon", ("code_length", [], None)),
        ("Metals", ("code_metallicity", ["metallicity"], None)),
        ("Phi", ("code_length", [], None)),
        ("FormationTime", ("code_time", ["creation_time"], None)),
    )



class TipsyFieldInfo(SPHFieldInfo):

    aux_particle_fields = {
            'uDotFB':("uDotFB", ("code_mass * code_velocity**2", ["uDotFB"], None)),
            'uDotAV':("uDotAV", ("code_mass * code_velocity**2", ["uDotAV"], None)),
            'uDotPdV':("uDotPdV", ("code_mass * code_velocity**2", ["uDotPdV"], None)),
            'uDotHydro':("uDotHydro", ("code_mass * code_velocity**2", ["uDotHydro"], None)),
            'uDotDiff':("uDotDiff", ("code_mass * code_velocity**2", ["uDotDiff"], None)),
            'uDot':("uDot", ("code_mass * code_velocity**2", ["uDot"], None)),
            'coolontime':("coolontime", ("code_time", ["coolontime"], None)),
            'timeform':("timeform", ("code_time", ["timeform"], None)),
            'massform':("massform", ("code_mass", ["massform"], None)),
            'HI':("HI", ("dimensionless", ["HI"], None)),
            'HII':("HII", ("dimensionless", ["HII"], None)),
            'HeI':("HeI", ("dimensionless", ["HeI"], None)),
            'HeII':("HeII", ("dimensionless", ["HeII"], None)),
            'OxMassFrac':("OxMassFrac", ("dimensionless", ["OxMassFrac"], None)),
            'FeMassFrac':("FeMassFrac", ("dimensionless", ["FeMassFrac"], None)),
            'c':("c", ("code_velocity", ["c"], None)),
            'acc':("acc", ("code_velocity / code_time", ["acc"], None)),
            'accg':("accg", ("code_velocity / code_time", ["accg"], None))}
    def __init__(self, pf, field_list, slice_info = None):
        for field in field_list:
            if field[1] in self.aux_particle_fields.keys() and \
                self.aux_particle_fields[field[1]] not in self.known_particle_fields:
                self.known_particle_fields += (self.aux_particle_fields[field[1]],)
        super(TipsyFieldInfo,self).__init__(pf, field_list, slice_info)


        

class OWLSFieldInfo(SPHFieldInfo):

    _species_fractions = ['H_fraction', 'He_fraction', 'C_fraction',
                          'N_fraction', 'O_fraction', 'Ne_fraction',
                          'Mg_fraction', 'Si_fraction', 'Fe_fraction']

    # override
    #--------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        
        new_particle_fields = (
            ('Hydrogen', ('', ['H_fraction'], None)),
            ('Helium', ('', ['He_fraction'], None)),
            ('Carbon', ('', ['C_fraction'], None)),
            ('Nitrogen', ('', ['N_fraction'], None)),
            ('Oxygen', ('', ['O_fraction'], None)),
            ('Neon', ('', ['Ne_fraction'], None)),
            ('Magnesium', ('', ['Mg_fraction'], None)),
            ('Silicon', ('', ['Si_fraction'], None)),
            ('Iron', ('', ['Fe_fraction'], None))
            )

        self.known_particle_fields += new_particle_fields
        
        super(OWLSFieldInfo,self).__init__( *args, **kwargs )


        
    def setup_fluid_fields(self):
        # here species_name is "H", "He", etc
        for s in self._species_fractions:
            species_name = s.split('_')[0]
            add_species_field_by_fraction(self, "gas", species_name)
