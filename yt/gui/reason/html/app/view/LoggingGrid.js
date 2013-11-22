/**********************************************************************
Logging Grid for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.LoggingGrid', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.logginggrid',
    title: 'Logging Output',
    store: 'LogEntries',
    itemId: 'logentries',
    columns: [ {header: 'Message', dataIndex:'record', flex:1} ],
    viewConfig: {
        cls: 'logentry',
    },
});
