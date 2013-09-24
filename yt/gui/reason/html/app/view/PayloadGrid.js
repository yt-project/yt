/**********************************************************************
Received Payloads Grid for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.PayloadGrid', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.payloadgrid',
    title: 'Received Payloads',
    store: 'Payloads',
    columns: [ {header: 'Payload Type', dataIndex:'payloadtype', flex: 1.0},
               {header: 'Widget Type', dataIndex: 'widgettype', flex: 1.0},
               {header: 'Widget Varname', dataIndex: 'varname', flex: 3.0},
    ],
});

