/**********************************************************************
The Phase Plot Creator

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

