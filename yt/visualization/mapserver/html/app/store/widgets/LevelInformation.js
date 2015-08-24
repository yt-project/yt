/**********************************************************************
Level Info store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.widgets.LevelInformation', {
    extend: 'Ext.data.Store',
    id: 'leveldata',
    fields: [
       {name: 'level', type:'int'},
       {name: 'cell_count', type: 'int'},
       {name: 'grid_count', type: 'int'},
       {name: 'grid_rel', type: 'float'},
       {name: 'cell_rel', type: 'float'},
    ],
    data: [],
});

