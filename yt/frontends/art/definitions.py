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
    'particle_type'
]

filename_pattern = {				
	'amr':'10MpcBox_csf512_%s.d',
	'particle_header':'PMcrd%s.DAT',
	'particle_data':'PMcrs0%s.DAT',
	'particle_stars':'stars_%s.dat'
}

amr_header_struct = [
    ('>i','pad byte'),
    ('>256s','jname'),
    ('>i','pad byte'),
    ('>i','pad byte'),
    ('>i','istep'),
    ('>d','t'),
    ('>d','dt'),
    ('>f','aexpn'),
    ('>f','ainit'),
    ('>i','pad byte'),
    ('>i','pad byte'),
    ('>f','boxh'),
    ('>f','Om0'),
    ('>f','Oml0'),
    ('>f','Omb0'),
    ('>f','hubble'),
    ('>i','pad byte'),
    ('>i','pad byte'),
    ('>i','nextras'),
    ('>i','pad byte'),
    ('>i','pad byte'),
    ('>f','extra1'),
    ('>f','extra2'),
    ('>i','pad byte'),
    ('>i','pad byte'),
    ('>256s','lextra'),
    ('>256s','lextra'),
    ('>i','pad byte'),
    ('>i', 'pad byte'),
    ('>i', 'min_level'),
    ('>i', 'max_level'),
    ('>i', 'pad byte'),
]

particle_header_struct =[
    ('>i','pad'),
    ('45s','header'), 
    ('>f','aexpn'),
    ('>f','aexp0'),
    ('>f','amplt'),
    ('>f','astep'),
    ('>i','istep'),
    ('>f','partw'),
    ('>f','tintg'),
    ('>f','Ekin'),
    ('>f','Ekin1'),
    ('>f','Ekin2'),
    ('>f','au0'),
    ('>f','aeu0'),
    ('>i','Nrow'),
    ('>i','Ngridc'),
    ('>i','Nspecies'),
    ('>i','Nseed'),
    ('>f','Om0'),
    ('>f','Oml0'),
    ('>f','hubble'),
    ('>f','Wp5'),
    ('>f','Ocurv'),
    ('>f','Omb0'),
    ('>%ds'%(396),'extras'),
    ('>f','unknown'),
    ('>i','pad')
]

constants = {
    "Y_p":0.245,
    "gamma":5./3.,
    "T_CMB0":2.726,
    "T_min":300.,
    "ng":128,
    "wmu":4.0/(8.0-5.0*0.245)
}

seek_extras = 137
