/**********************************************************************
Widget Types Grid for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.WidgetTypesGrid', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.widgettypesgrid',
    title: 'Known Widget Types',
    store: 'WidgetTypes',
    columns: [ {header: 'Short Name', dataIndex:'widgetname'},
               {header: 'Display Name', dataIndex: 'displayname'},
               {header: 'PFs', dataIndex: 'pfs'},
               {header: 'DOs', dataIndex: 'objs'}
    ],
});

