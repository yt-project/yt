/**********************************************************************
Logging entry controller for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.controller.Logging', {
    extend: 'Ext.app.Controller',
    stores: [ 'LogEntries' ],
    view: ['LoggingGrid'],

    refs: [
        { ref: 'logEntries',
          selector: '#logentries'
        },
    ],

    init: function() {
        this.application.addListener({
            payloadlogentry: {fn: this.addLogPayload, scope: this},
            logentry: {fn: this.addLogEntry, scope: this},
        });
    },

    addLogPayload: function(payload) {
        this.addLogEntry(payload['log_entry']);
    },

    addLogEntry: function(text) {
        this.getLogEntriesStore().add({record: text});
        var i = this.getLogEntriesStore().getCount();
        this.getLogEntries().getView().focusRow(i - 1);
    },
});
