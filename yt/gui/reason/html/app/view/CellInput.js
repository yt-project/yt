/**********************************************************************
Cell input in Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.CellInput', {
    extend: 'Ext.panel.Panel',
    requires: ['Ext.form.Panel'],
    alias: 'widget.notebookinput',
    title: 'yt Input',
    width: '100%',
    flex: 0.3,
    layout: {type: 'hbox',
             pack: 'start',
             align: 'stretch',
             },
    items: [{
                xtype: 'form',
                flex: 1,
                width: '100%',
                layout: 'fit',
                items: [
                    {
                      xtype: 'textarea',
                      id: 'inputline',
                      autoScroll: true,
                      name: 'line',
                      width:'100%',
                      allowBlank: 'True',
                      fieldCls: 'inputLine',
                    },
                ],
            }, {
                xtype: 'button',
                iconCls: 'doubledownarrow',
                tooltip: 'Execute Cell',
                id: 'executecellbutton',
                width: 24,
                height: '100%',
            }],
});
