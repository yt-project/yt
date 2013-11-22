/**********************************************************************
The Isocontour Creator

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.IsocontourCreator", {
    extend: 'Ext.window.Window',
    alias: 'widget.isocontourcreator',
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
                fieldLabel: 'Field',
                id: 'field',
                width: 290,
                queryMode: 'local',
                editable: false,
                triggerAction: 'all',
                value: 'Density'
              },{
                xtype:'checkboxfield',
                fieldLabel: 'Relative Value?',
                itemId: 'relValue',
                value: true,
                checked: true,
                width: 290,
                allowBlank:false,
              },{
                xtype:'textfield',
                fieldLabel: 'Value',
                itemId: 'value',
                value: '0.5',
                width: 290,
                allowBlank:false,
              }
            ],
            buttons: [
                {
                    text: 'Extract',
                    itemId: 'extract',
                },{
                    text: 'Cancel',
                    itemId: 'cancel',
                }
            ]
        }],
});

