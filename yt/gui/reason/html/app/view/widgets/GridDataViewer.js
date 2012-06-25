/**********************************************************************
The Grid Data Viewer View

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

Ext.define("Reason.view.widgets.GridDataViewer", {
    extend: 'Ext.grid.Panel',
    title: 'This should not be visible.',
    alias: 'widget.griddataviewer',
    iconCls: 'graph',
    autoScroll: true,
    width: '100%',
    height: '100%',
    closable: true,
    columns: [
         {
             id: 'grid_id',
             header: 'Grid ID',
             width: 100,
             dataIndex: 'grid_id',
             sortable: true,
         }, {
             id: 'left_edge_x',
             header: 'Left Edge x',
             width: 100,
             dataIndex: 'left_edge_x',
             sortable: true,
         }, {
             id: 'left_edge_y',
             header: 'Left Edge y',
             width: 100,
             dataIndex: 'left_edge_y',
             sortable: true,
         }, {
             id: 'left_edge_z',
             header: 'Left Edge z',
             width: 100,
             dataIndex: 'left_edge_z',
             sortable: true,
         }, {
             id: 'right_edge_x',
             header: 'Right Edge x',
             width: 100,
             dataIndex: 'right_edge_x',
             sortable: true,
         }, {
             id: 'right_edge_y',
             header: 'Right Edge y',
             width: 100,
             dataIndex: 'right_edge_y',
             sortable: true,
         }, {
             id: 'right_edge_z',
             header: 'Right Edge z',
             width: 100,
             dataIndex: 'right_edge_z',
             sortable: true,
         }, {
             id: 'dim_x',
             header: 'DimX',
             width: 100,
             dataIndex: 'dim_x',
             sortable: true,
         }, {
             id: 'dim_y',
             header: 'DimY',
             width: 100,
             dataIndex: 'dim_y',
             sortable: true,
         }, {
             id: 'dim_z',
             header: 'DimZ',
             width: 100,
             dataIndex: 'dim_z',
             sortable: true,
         }, {
             id: 'cells',
             header: 'Cells',
             width: 100,
             dataIndex: 'cells',
             sortable: true,
         },
     ],
});

