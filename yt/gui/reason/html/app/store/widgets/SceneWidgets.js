/**********************************************************************
Widget Store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.widgets.SceneWidgets', {
    extend: 'Ext.data.Store',
    id: 'scenewidgets',
    fields: [
       {name: 'name', type: 'string'},
       {name: 'enabled', type: 'boolean', defaultValue: true},
       {name: 'type', type: 'string'},
       {name: 'widget'},
    ],
    data: [],
});
