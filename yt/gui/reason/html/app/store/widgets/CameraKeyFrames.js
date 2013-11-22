/**********************************************************************
Camera Key Frame store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.widgets.CameraKeyFrames', {
    extend: 'Ext.data.Store',
    id: 'camerakeyframes',
    fields: [
       {name: 'time', type:'float'},
       {name: 'pos_x', type:'float'},
       {name: 'pos_y', type:'float'},
       {name: 'pos_z', type:'float'},
       {name: 'view'},
    ],
    data: [],
});

