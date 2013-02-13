"""
Definitions specific to ART

Author: Christopher E. Moody <cemoody@ucsc.edu>
Affiliation: UC Santa Cruz
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Christopher E. Moody.  All Rights
  Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

fluid_fields= [ 
    'Density',
    'TotalEnergy',
    'XMomentumDensity',
    'YMomentumDensity',
    'ZMomentumDensity',
    'Pressure',
    'Gamma',
    'GasEnergy',
    'MetalDensitySNII',
    'MetalDensitySNIa',
    'PotentialNew',
    'PotentialOld'
]

hydro_struct = [('pad1',1,'>i'),('idc',1,'>i'),('iOctCh',1,'>i')]
for field in fluid_fields:
    hydro_struct += (field,1,'>f'),
hydro_struct += ('pad2',1,'>i'),

particle_fields= [
    'particle_age',
    'particle_index',
    'particle_mass',
    'particle_mass_initial',
    'particle_creation_time',
    'particle_metallicity1',
    'particle_metallicity2',
    'particle_metallicity',
    'particle_position_x',
    'particle_position_y',
    'particle_position_z',
    'particle_velocity_x',
    'particle_velocity_y',
    'particle_velocity_z',
    'particle_type',
    'particle_index'
]

particle_star_fields = [
    'particle_age',
    'particle_mass',
    'particle_mass_initial',
    'particle_creation_time',
    'particle_metallicity1',
    'particle_metallicity2',
    'particle_metallicity',
]

filename_pattern = {				
	'amr':'10MpcBox_csf512_%s.d',
	'particle_header':'PMcrd%s.DAT',
	'particle_data':'PMcrs0%s.DAT',
	'particle_stars':'stars_%s.dat'
}

filename_pattern_hf = {				
	'particle_header':'PMcrd_%s.DAT',
	'particle_data':'PMcrs0_%s.DAT',
}

amr_header_struct = [
    ('pad byte',1,'>i'),
    ('jname',1,'>256s'),
    ('pad byte',1,'>i'),
    ('pad byte',1,'>i'),
    ('istep',1,'>i'),
    ('t',1,'>d'),
    ('dt',1,'>d'),
    ('aexpn',1,'>f'),
    ('ainit',1,'>f'),
    ('pad byte',1,'>i'),
    ('pad byte',1,'>i'),
    ('boxh',1,'>f'),
    ('Om0',1,'>f'),
    ('Oml0',1,'>f'),
    ('Omb0',1,'>f'),
    ('hubble',1,'>f'),
    ('pad byte',1,'>i'),
    ('pad byte',1,'>i'),
    ('nextras',1,'>i'),
    ('pad byte',1,'>i'),
    ('pad byte',1,'>i'),
    ('extra1',1,'>f'),
    ('extra2',1,'>f'),
    ('pad byte',1,'>i'),
    ('pad byte',1,'>i'),
    ('lextra',1,'>256s'),
    ('lextra',1,'>256s'),
    ('pad byte',1,'>i'),
    ( 'pad byte',1,'>i'),
    ( 'min_level',1,'>i'),
    ( 'max_level',1,'>i'),
    ( 'pad byte',1,'>i')
]

particle_header_struct =[
    ('pad',1,'>i'),
    ('header',1,'45s', 
    ('aexpn',1,'>f'),
    ('aexp0',1,'>f'),
    ('amplt',1,'>f'),
    ('astep',1,'>f'),
    ('istep',1,'>i'),
    ('partw',1,'>f'),
    ('tintg',1,'>f'),
    ('Ekin',1,'>f'),
    ('Ekin1',1,'>f'),
    ('Ekin2',1,'>f'),
    ('au0',1,'>f'),
    ('aeu0',1,'>f'),
    ('Nrow',1,'>i'),
    ('Ngridc',1,'>i'),
    ('Nspecies',1,'>i'),
    ('Nseed',1,'>i'),
    ('Om0',1,'>f'),
    ('Oml0',1,'>f'),
    ('hubble',1,'>f'),
    ('Wp5',1,'>f'),
    ('Ocurv',1,'>f'),
    ('Omb0',1,'>f'),
    ('extras',1,'>%ds'%(396)),
    ('unknown',1,'>f'),
    ('pad',1,'>i')
]

star_struct = [
        ('>d',('tdum','adum')),
        ('>i','nstars'),
        ('>d',('ws_old','ws_oldi')),
        ('>f','mass'),
        ('>f','imass'),
        ('>f','tbirth'),
        ('>f','metallicity1'),
        ('>f','metallicity2')
        ]

star_name_map = {
        'particle_mass':'mass',
        'particle_mass_initial':'imass',
        'particle_age':'tbirth',
        'particle_metallicity1':'metallicity1',
        'particle_metallicity2':'metallicity2',
        'particle_metallicity':'metallicity',
        }

constants = {
    "Y_p":0.245,
    "gamma":5./3.,
    "T_CMB0":2.726,
    "T_min":300.,
    "ng":128,
    "wmu":4.0/(8.0-5.0*0.245)
}

seek_extras = 137
