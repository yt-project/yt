/**********************************************************************
File Listing store for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.FileListing', {
    extend: 'Ext.data.ArrayStore',
    id: 'filelisting',
    fields: [
        'filename', 
        {name: 'size', type: 'float'},
        'type'
    ],
    data: [],
});

