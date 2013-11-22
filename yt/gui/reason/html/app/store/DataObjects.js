/**********************************************************************
Logging entry store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.DataObjects', {
    extend: 'Ext.data.TreeStore',
    requires: ['Reason.model.DataObject'],
    model: 'Reason.model.DataObject',
    autoLoad: false,
    proxy: {
        type: 'memory',
        reader: {
            type: 'array',
        }
    },
});

