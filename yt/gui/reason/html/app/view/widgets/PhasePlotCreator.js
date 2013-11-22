/**********************************************************************
The Phase Plot Creator

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.PhasePlotCreator", {
    extend: 'Ext.window.Window',
    alias: 'widget.phaseplotcreator',
    width: 370,
    height: 200,
    modal: true,
    resizable: false,
    draggable: false,
    border: false,
    layout: 'fit',
    title: "Create Phase Plot",

    items: [{ xtype: 'form', // FormPanel
              labelWidth:80,
              frame:true,
              items: [ {
                xtype:'combo',
                fieldLabel: 'X Field',
                id: 'x_field',
                width: 350,
                queryMode: 'local',
                editable: false,
                triggerAction: 'all',
                value: 'Density'
            },{
                xtype:'combo',
                fieldLabel: 'Y Field',
                id: 'y_field',
                width: 350,
                queryMode: 'local',
                editable: false,
                triggerAction: 'all',
                value: 'Density'
            },{
                xtype:'combo',
                fieldLabel: 'Z Field',
                id: 'z_field',
                width: 350,
                queryMode: 'local',
                editable: false,
                triggerAction: 'all',
                value: 'Density'
            },{
                xtype:'combo',
                fieldLabel: 'Weight by',
                id: 'weight',
                width: 350,
                allowBlank:false,
                triggerAction: 'all',
                value: 'Just Sum',
                queryMode: 'local',
                store: ['Just Sum', 'Mass', 'Volume'],
            }],
            buttons: [
                {
                    text: 'Calculate',
                    itemId: 'create',
                },{
                    text: 'Cancel',
                    itemId: 'cancel',
                }
            ]
        }],
});

