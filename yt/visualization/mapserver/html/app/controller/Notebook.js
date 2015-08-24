/**********************************************************************
Notebook controller for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
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

    addRequest: function(request_id, command) {
        /*console.log("Adding request " + request_id);*/
        this.getRequestsStore().add({
            request_id: request_id, command: command,
        });
        reason.pending.update([this.getRequestsStore().count()]);
    },

    addCell: function(cell) {
        this.application.fireEvent("wipeinput");
        this.application.fireEvent("allowinput");
        if (cell['result_id'] != null) {
            var ind = this.getRequestsStore().find(
                'request_id', cell['result_id']);
            if (ind != -1) {
                /*console.log("Removing request " + cell['result_id']);*/
                var rec = this.getRequestsStore().removeAt(ind);
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
        this.getInputLine().setValue(payload['value']);
        this.getInputLine()._lock = true;
    },

    executeCell: function(line) {
        this.application.fireEvent("blockinput");
        reason.server.execute(line, false);
    },

    scrollToBottom: function() {
        var i = this.getCellValuesStore().getCount();
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

