/**********************************************************************
Logging entry store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.CellValues', {
    extend: 'Ext.data.Store',
    id: 'cellvalues',
    fields: ['input', 'output', 'raw_input', 'executiontime', 
        { name: 'image_data', type: 'string', defaultValue: '' },
        'result_id',
    ],
    data: [],
});

