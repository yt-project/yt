/**********************************************************************
Data Object Tree for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.DataObjectTree', {
    extend: 'Ext.tree.Panel',
    alias: 'widget.dataobjecttree',
    store: 'DataObjects',
    itemId: 'dataobjects',
    rootVisible: false,
    iconCls: 'nav',
    displayField: 'name',
    header: false,
    root: {
        expanded: true,
        text: '',
        leaf: false,
    },
    columns: {
        items: [ 
            { xtype: 'treecolumn',
              text: 'Object',
              sortable: false,
              dataIndex: 'name',
              flex: 1.0,
            },
        ]
    },
});

