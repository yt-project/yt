/**********************************************************************
The Plot Window Widget

Author: Cameron Hummels <chummels@gmail.com>
Affiliation: Columbia
Author: Jeffrey S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
                  store: ['None'],
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
                      text: 'Project', itemId: 'create',
                  },{
                      text: 'Cancel', itemId: 'cancel',
                  }
              ]
    }],
});

