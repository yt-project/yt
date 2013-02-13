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

#If not otherwise specified, we are big endian
endian = '>'

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
    ('jname',1,'256s'),
    (('istep','t','dt','aexpn','ainit'),1,'iddff'),
    (('boxh','Om0','Oml0','Omb0','hubble'),5,'f'),
    ('nextras',1,'i'),
    (('extra1','extra2'),2,'f'),
    ('lextra',1,'512s'),
    (('min_level','max_level'),2,'i')
]

particle_header_struct =[
    (('header',
     'aexpn','aexp0','amplt','astep',
     'istep',
     'partw','tintg',
     'Ekin','Ekin1','Ekin2',
     'au0','aeu0',
     'Nrow','Ngridc','Nspecies','Nseed',
     'Om0','Oml0','hubble','Wp5','Ocurv','Omb0',
     'extras','unknown'),
      1,
     '45sffffi'+'fffffff'+'iiii'+'ffffff'+'396s'+'f')
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
