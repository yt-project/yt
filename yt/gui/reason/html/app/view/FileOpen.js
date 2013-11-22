/**********************************************************************
File Open Dialog for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
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


