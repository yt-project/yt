/**********************************************************************
Notebook controller for Reason

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

Ext.define('Reason.controller.Notebook', {
    extend: 'Ext.app.Controller',
    stores: [ 'CellValues' , 'Requests'],
    views: ['Notebook', 'CellView'],
    refs: [
        { ref: 'inputLine',
          selector: '#inputline'
        },
        { ref: 'cellDisplay',
          selector: 'notebookcells#cells',
          xtype: 'notebookcells',
          autoCreate: true,
          itemId: 'cells',
        },
    ],

    init: function() {
        this.application.addListener({
            payloadcell: {fn: this.addCell, scope: this},
            payloadscript: {fn: this.loadScript, scope: this},
            executecell: {fn: this.executeCell, scope: this},
            wipeinput:   {fn: this.wipeInputLine, scope: this},
            blockinput:  {fn: this.blockInput, scope: this},
            allowinput:  {fn: this.allowInput, scope: this},
            scrolltobottom: {fn: this.scrollToBottom, scope: this},
        })
        this.control({
            '#executecellbutton': {
                click: function(f, e) {
                    this.executeCell(this.getInputLine().getValue());
                }
            },
            '#inputline': {
                specialkey: function(field, e, opts){
                    if (e.shiftKey && e.getKey() == e.ENTER) {
                        this.executeCell(field.getValue());
                    }
                },
            },
        });
        this.callParent(arguments);
    },

    addRequest: function(request_id) {
        this.getRequestsStore().add({
            request_id: request_id, pending: true,
        });
    },

    addCell: function(cell) {
        this.application.fireEvent("wipeinput");
        this.application.fireEvent("allowinput");
        if (cell['result_id'] != null) {
            var ind = this.getRequestsStore().find(
                'request_id', cell['result_id']);
            if (ind != -1) {
                var rec = this.getRequestsStore().getAt(ind);
                rec.data.pending = false;
            }
            reason.pending.update([this.getRequestsStore().count()]);
        }
        if (cell['hide'] == true) { return; }
        this.getCellValuesStore().add({
            input: cell['input'],
            output: cell['output'],
            raw_input: cell['raw_input'],
            image_data: cell['image_data'],
            executiontime: cell['executiontime'],
            result_id: cell['result_id'],
        });
        this.application.fireEvent("scrolltobottom");
    },

    loadScript: function(payload) {
        console.log("Loading script ...");
        this.getInputLine().setValue(payload['value']);
        this.getInputLine()._lock = true;
    },

    executeCell: function(line) {
        this.application.fireEvent("blockinput");
        console.log("Asked to execute " + line);
        reason.server.execute(line, false);
    },

    scrollToBottom: function() {
        var i = this.getCellValuesStore().getCount();
        console.log("Scrolling to bottom: " + i);
        this.getCellDisplay().getView().focusRow(i-1);
    },
    
    wipeInputLine: function() {
        if(this.getInputLine()._lock == true) {
            this.getInputLine()._lock = false;
        } else {
            this.getInputLine().setValue("");
        }
    },

    blockInput: function() {
        this.getInputLine().addClass("cell_waiting");
        this.getInputLine().setReadOnly(true);
    },

    allowInput: function() {
        this.getInputLine().removeCls("cell_waiting");
        this.getInputLine().setReadOnly(false);
        var application = this.application;
    },

});

