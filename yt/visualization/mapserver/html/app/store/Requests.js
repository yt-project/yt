/**********************************************************************
A store for outstanding requests

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.store.Requests', {
    extend: 'Ext.data.Store',
    id: 'requestsstore',
    fields: ['request_id', 'command'],
    data: [],
});

