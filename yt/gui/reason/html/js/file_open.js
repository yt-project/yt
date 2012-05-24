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

    function fillStore(f, a){
        if(a.status == false){
          Ext.Msg.alert("Error", "Something has gone wrong.");
          examine = {f: f, a: a};
          return;
        }
        examine = a;
        filestore.removeAll();
        var rec = [];
        filestore.loadData(a.result['objs']);
    }

    var win = new Ext.Window({
        layout:'fit',
        width:480,
        height:320,
        modal:true,
        resizable:true,
        draggable:true,
        border:false,
        title:'Open File',
        items: [
            { xtype:'listview',
              store: filestore ,
              emptyText: 'No images to display',
              columns: [
              {
                  header: 'Type',
                  width: 0.1,
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
              }]
            },
        ]
    });
    yt_rpc.ExtDirectREPL.file_listing(
          {base_dir:"", sub_dir:""}, fillStore);
    win.show(this);
}
