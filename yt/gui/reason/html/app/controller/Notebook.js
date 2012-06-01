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
    stores: [ 'CellValues' ],
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

    addCell: function(cell) {
        this.application.fireEvent("wipeinput");
        this.application.fireEvent("allowinput");
        this.getCellValuesStore().add({
            input: cell['input'],
            output: cell['output'],
            raw_input: cell['raw_input'],
            executiontime: cell['executiontime'],
        });
        this.application.fireEvent("scrolltobottom");
    },
    executeCell: function(line) {
        this.application.fireEvent("blockinput");
        console.log("Asked to execute " + line);
        yt_rpc.ExtDirectREPL.execute({code:line}, this.cellExecuted);
    },

    scrollToBottom: function() {
        var i = this.getCellValuesStore().getCount();
        console.log("Scrolling to bottom: " + i);
        this.getCellDisplay().getView().focusRow(i-1);
    },
    
    wipeInputLine: function() {
        this.getInputLine().setValue("");
    },

    blockInput: function() {
        this.getInputLine().addClass("cell_waiting");
        this.getInputLine().setReadOnly(true);
    },

    allowInput: function() {
        this.getInputLine().removeCls("cell_waiting");
        this.getInputLine().setReadOnly(false);
        console.log("Calling FileList");
        var application = this.application;
    },

    cellExecuted: function(result) {
        console.log("Cell Executed!");
    },
});

