/**********************************************************************
A file opener

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

function open_file() {
    var filestore = new Ext.data.ArrayStore({
      fields: ['filename', 
               {name:'size', type:'float'},
               'type'
      ]
    });
    var cur_dir;
    function fillStore(f, a){
        if(a.status == false){
          Ext.Msg.alert("Error", "Something has gone wrong.");
          return;
        }
        if(a.result['change'] == false) {
          win.get("current_file").setValue(cur_dir);
          return;
        }
        filestore.removeAll();
        var rec = [];
        filestore.loadData(a.result['objs']);
        cur_dir = a.result['cur_dir'];
        win.get("current_file").setValue(cur_dir);
    }

    var win = new Ext.Window({
        layout:'vbox',
        layoutConfig: {
            align: 'stretch',
            pack: 'start',
            defaultMargins: "5px 5px 5px 5px",
        },
        width:540,
        height:480,
        modal:true,
        resizable:true,
        draggable:true,
        title:'Open File',
        items: [
            { xtype: 'textfield',
              id: 'current_file',
              listeners: {
                specialkey: function(f, e) {
                  if (e.getKey() != e.ENTER) { return; }
                  yt_rpc.ExtDirectREPL.file_listing(
                        {base_dir:f.getValue(), sub_dir:''}, fillStore);
                }
              }
            }, {
              xtype:'listview',
              id: 'file_listing',
              store: filestore ,
              singleSelect:true,
              emptyText: 'No images to display',
              flex: 1.0,
              columns: [
              {
                  header: 'Type',
                  width: 0.1,
                  tpl: '<img src="images/file_dialog_{type}.png" width=16 height=16>',
                  dataIndex: 'type'
              },{
                  header: 'Filename',
                  width: .75,
                  dataIndex: 'filename'
              },{
                  header: 'Size',
                  dataIndex: 'size',
                  tpl: '{size:fileSize}',
                  align: 'right',
                  cls: 'listview-filesize'
              }],
              listeners: {
                dblclick: function(view, index, node, e) {
                    var fileRecord = filestore.getAt(index).data;
                    if (fileRecord.type == 'directory') {
                      yt_rpc.ExtDirectREPL.file_listing(
                            {base_dir:cur_dir, sub_dir:fileRecord.filename},
                            fillStore);
                    } else {
                      yt_rpc.ExtDirectREPL.load(
                            {base_dir:cur_dir, filename:fileRecord.filename},
                            handle_result);
                      win.destroy();
                    }
                },
                selectionchange: function(view, index, node, e) {
                },
              },
            }, {
              xtype: 'panel',
              height: 40,
              layout: 'hbox',
              layoutConfig: {
                  align: 'stretch',
                  pack: 'start',
                  defaultMargins: "5px 5px 5px 5px",
              },
              items: [
                { flex: 1.0, xtype: 'button', text: 'Cancel',
                    handler: function(b, e) { win.destroy(); } },
                { flex: 1.0, xtype: 'button', text: 'Load',
                    handler: function(b, e) {
                      filename = "";
                      var fl = win.get("file_listing");
                      if (fl.getSelectionCount() == 1) {
                        filename = fl.getSelectedRecords()[0].data.filename;
                      }
                      yt_rpc.ExtDirectREPL.load(
                            {base_dir:cur_dir, filename:filename},
                            handle_result);
                      win.destroy();
                    }
                },
              ],
            },
        ],
    });
    yt_rpc.ExtDirectREPL.file_listing(
          {base_dir:"", sub_dir:""}, fillStore);
    win.show(this);
}
