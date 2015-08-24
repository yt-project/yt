/**********************************************************************
Grid data store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
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

