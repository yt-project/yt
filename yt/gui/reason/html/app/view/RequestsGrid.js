/**********************************************************************
Requests Grid for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.RequestsGrid', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.requestsgrid',
    title: 'Pending Requests',
    store: 'Requests',
    columns: [ {header: 'Request ID', dataIndex:'result_id', width:80},
               {header: 'Command', dataIndex: 'command', flex: 1.0},
    ],
});
