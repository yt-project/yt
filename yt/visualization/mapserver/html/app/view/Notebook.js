/**********************************************************************
Notebook view for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.Notebook', {
    extend: 'Ext.panel.Panel',
    requires: [ 'Reason.view.CellView',
                'Reason.view.CellInput' ],
    alias: 'widget.ytnotebook',
    title: 'yt Notebook',
    width: '100%',
    layout: {
        type: 'vbox',
        align:'stretch',
        pack:'start',
    },
    closable: false,
    autoScroll: false,
    iconCls: 'console',
    items: [
        { xtype: 'notebookinput', },
        { xtype: 'splitter' },
        { xtype: 'notebookcells', },
    ],
});


