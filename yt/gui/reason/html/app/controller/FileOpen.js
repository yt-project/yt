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

Ext.define('Reason.controller.FileOpen', {
    extend: 'Ext.app.Controller',
    stores: ['FileListing'],
    requires: ["Reason.view.FileOpen"],
    
    init: function() {

    },

    openDialog: function() {
        this.currentDirectory = '',
        this.getFileListingStore().removeAll();
        reason.server.method("file_listing",
            {base_dir:"", sub_dir:""}, this.fillStore);
        this.win = Ext.widget("fileopenwindow");
        this.win.query("#file_listing")[0].bindStore(this.getFileListingStore());
        /* Set up our handlers */
        this.control({
            '#current_file': { specialkey: this.specialKeyPressed, },
            '#file_listing': { celldblclick: this.doubleClicked, },
            '#load': { click: this.loadFile, },
            '#cancel': { click: this.cancelWindow, }
        });
        this.win.show();
    },

    specialKeyPressed: function(f, e) {
        if (e.getKey() != e.ENTER) { return; }
        reason.server.method("file_listing",
            {base_dir:f.getValue(), sub_dir:""}, this.fillStore);
    },  

    doubleClicked: function(view, td, cellIndex, record) {
        conf = {'base_dir': this.currentDirectory};
        if (record.data.type == 'directory') {
          conf['sub_dir'] = record.data.filename;
          reason.server.method('file_listing', conf, this.fillStore);
        } else {
          conf['filename'] = record.data.filename;
          reason.server.method('load', conf);
          this.win.destroy();
        }
    },

    loadFile: function(b, e) {
        var conf = {filename: '', base_dir: this.currentDirectory};
        var fl = this.win.query("#file_listing")[0];
        var idx = fl.getSelectionModel().getSelection();
        if (idx == 1) {
            conf['filename'] = this.getFileListingStore().getAt(idx).data.filename;
        }
        reason.server.method('load', conf);
    },

    cancelWindow:  function(b, e) {
        this.win.destroy();
    },

    fillStore: function(f, a) {
        var con = reason.getController("FileOpen");
        if(a.status == false){
          Ext.Msg.alert("Error", "Something has gone wrong.");
          con.window.destroy();
        }
        if(a.result['change'] == false) {
          con.win.query("#current_file")[0].setValue(con.currentDirectory);
          return;
        }
        con.getFileListingStore().removeAll();
        var rec = [];
        con.getFileListingStore().loadData(a.result['objs']);
        con.currentDirectory = a.result['cur_dir'];
        con.win.query("#current_file")[0].setValue(con.currentDirectory);
    },
});
