/**********************************************************************
The Isocontour Creator

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

