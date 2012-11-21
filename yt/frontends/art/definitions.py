"""
Definitions specific to ART

Author: Christopher E. Moody <cemoody@ucsc.ed>
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

art_particle_field_names = [
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
'particle_type']

filename_pattern = {				
	'amr':'10MpcBox_csf512_%s.d',
	'particle_header':'PMcrd%s.DAT',
	'particle_data':'PMcrs0%s.DAT',
	'particle_stars':'stars_%s.dat'
}
