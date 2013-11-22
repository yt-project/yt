/**********************************************************************
The Plot Window Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.PlotWindowCreator", {
    extend: 'Ext.window.Window',
    alias: 'widget.plotwindowcreator',
    width: 480,
    height: 200,
    modal: true,
    resizable: true,
    draggable: false,
    border: false,
    layout: 'fit',
    title: "You shouldn't see this.",

    items: [{ xtype: 'form',
              labelWidth:80,
              frame:true,
              defaults: {
                  bodyStyle:'padding:15px 20px'
              },
              layout: {
                  type: 'table',
                  columns: 2
              },
              items: [{
                  xtype:'combo',
                  fieldLabel: 'Type of Plot',
                  itemId: 'plotType',
                  store: ['Slice', 'Projection'],
                  width: 200,
                  allowBlank: false,
                  editable: false,
                  value: 'Slice'
              },{
                  xtype:'checkbox',
                  fieldLabel: 'Center on Max',
                  itemId: 'maxDens',
                  width: 200,
                  allowBlank:false,
              },{
                  xtype:'combo',
                  fieldLabel: 'Axis',
                  itemId: 'axis',
                  store:['X','Y','Z'],
                  width: 200,
                  allowBlank: false,
                  editable: false,
                  triggerAction: 'all',
                  value: 'X',
              },{
                  xtype:'textfield',
                  fieldLabel: 'Center X',
                  itemId: 'slice_x_center',
                  value: '0.5',
                  width: 200,
                  allowBlank:false,
              },{
                  xtype:'combo',
                  fieldLabel: 'Field',
                  itemId: 'field',
                  store: [],
                  displayField: 'text',
                  valueield: 'text',
                  width: 200,
                  allowBlank:false,
                  triggerAction: 'all',
                  value: 'Density'
              },{
                  xtype:'textfield',
                  fieldLabel: 'Center Y',
                  itemId: 'slice_y_center',
                  value: '0.5',
                  width: 200,
                  allowBlank:false,
              },{
                  xtype:'combo',
                  fieldLabel: 'Weight Field',
                  itemId: 'weightField',
                  store: [{text:'None'},],
                  displayField: 'text',
                  valueield: 'text',
                  width: 200,
                  allowBlank:false,
                  triggerAction: 'all',
                  value: 'None'
              },{
                  xtype:'textfield',
                  fieldLabel: 'Center Z',
                  itemId: 'slice_z_center',
                  value: '0.5',
                  width: 200,
                  allowBlank:false,
              }],
              buttons: [
                  {
                      text: 'Go', itemId: 'create',
                  },{
                      text: 'Cancel', itemId: 'cancel',
                  }
              ]
    }],
});

