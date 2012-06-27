/**********************************************************************
File Open Dialog for Reason

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

Ext.define('Reason.view.FileOpen', {
    extend: 'Ext.window.Window',
    alias: 'widget.fileopenwindow',

    layout: {
        type: 'vbox',
        align: 'stretch',
        pack: 'start',
        defaultMargins: '5px 5px 5px 5px',
    },
    width: 540,
    height: 480,
    modal:true,
    resizable:true,
    draggable:true,
    title:'Open File',

    items: [
        {
          xtype: 'textfield',
          itemId: 'current_file',
          height: 35,
        }, {
          xtype:'gridpanel',
          itemId: 'file_listing',
          singleSelect:true,
          emptyText: 'No files to display',
          flex: 1.0,
          columns: [
          {
              header: 'Type',
              width: 24,
              renderer: function(value) {
                return Ext.String.format('<img src="reason/resources/images/' + 
                                         'file_dialog_{0}.png"' + 
                                         'width=16 height=16>', value);
              },
              dataIndex: 'type'
          },{
              header: 'Filename',
              flex: 3.0,
              dataIndex: 'filename'
          },{
              header: 'Size',
              dataIndex: 'size',
              tpl: '{size:fileSize}',
              align: 'right',
              cls: 'listview-filesize',
              flex: 1.0,
          }],
        }, {
          xtype: 'panel',
          height: 40,
          layout: {
              type: 'hbox',
              align: 'stretch',
              pack: 'start',
              defaultMargins: "5px 5px 5px 5px",
          },
          items: [
            { flex: 1.0,
              xtype: 'button',
              text: 'Cancel', 
              itemId: 'cancel',
            }, { 
              flex: 1.0,
              xtype: 'button',
              text: 'Load',
              itemId: 'load',
            },
          ],
        },
    ],
});


