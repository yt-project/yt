/**********************************************************************
Camera Path Element store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.widgets.CameraPathElements', {
    extend: 'Ext.data.Store',
    id: 'camerapathelements',
    fields: [
       {name: 'position'},
       {name: 'focus'},
       {name: 'up'},
    ],
    data: [],
});

