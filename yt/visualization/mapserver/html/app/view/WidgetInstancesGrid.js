/**********************************************************************
Widget Instances Grid for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.WidgetInstancesGrid', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.widgetinstancesgrid',
    title: 'Existing Widgets',
    store: 'WidgetInstances',
    columns: [ {header: 'Widget Type', dataIndex:'widgettype', flex: 1.0},
               {header: 'Widget ID', dataIndex: 'widgetid', flex: 3.0},
    ],
});

