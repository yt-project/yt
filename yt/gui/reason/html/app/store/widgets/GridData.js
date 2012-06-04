/**********************************************************************
Grid data store for Reason

Author: Cameron Hummels <chummels@gmail.com>
Affiliation: Columbia
Author: Jeffrey S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
***********************************************************************/

Ext.define('Reason.store.widgets.GridData', {
    extend: 'Ext.data.Store',
    id: 'griddata',
    fields: [
       {name: 'grid_id', type:'int'},
       {name: 'level', type:'int'},
       {name: 'left_edge_x', type: 'float'},
       {name: 'left_edge_y', type: 'float'},
       {name: 'left_edge_z', type: 'float'},
       {name: 'right_edge_x', type: 'float'},
       {name: 'right_edge_y', type: 'float'},
       {name: 'right_edge_z', type: 'float'},
       {name: 'dim_x', type: 'int'},
       {name: 'dim_y', type: 'int'},
       {name: 'dim_z', type: 'int'},
       {name: 'cells', type: 'int'},
    ],
    data: [],
});

