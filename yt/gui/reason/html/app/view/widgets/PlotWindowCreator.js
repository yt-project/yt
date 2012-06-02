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
    layout: 'fit',
    width: 370,
    height: 220,
    modal: true,
    resizable: false,
    draggable: false,
    border: false,
    title: "You shouldn't see this.",

    items: [{ xtype: 'form',
              labelWidth:80,
              frame:true,
              items: [{
                  xtype:'combo',
                  fieldLabel: 'Axis',
                  itemId: 'axis',
                  store:['X','Y','Z'],
                  width: 230,
                  allowBlank:false,
                  triggerAction: 'all',
                  value: 'X',
              },{
                  xtype:'checkbox',
                  fieldLabel: 'Center on Max',
                  itemId: 'maxDens',
                  width: 90,
                  allowBlank:false,
              },{
                  xtype:'combo',
                  fieldLabel: 'Field',
                  itemId: 'field',
                  store: [],
                  width: 230,
                  allowBlank:false,
                  triggerAction: 'all',
                  value: 'Density'
              },{
                  xtype:'combo',
                  fieldLabel: 'Weight Field',
                  itemId: 'weightField',
                  store: ['None'],
                  width: 230,
                  allowBlank:false,
                  triggerAction: 'all',
                  value: 'None'
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

