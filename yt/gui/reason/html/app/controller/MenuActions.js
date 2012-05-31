/**********************************************************************
Menu actions in Reason

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


Ext.define('Reason.controller.MenuActions', {
    extend: 'Ext.app.Controller',
    view: ['MainMenu'],

    init: function() {
        this.callParent(arguments);
    },

    addLogEntry: function(payload) {
        examine = payload;
        this.getLogEntriesStore().add(
            {record: payload['log_entry']}
        );
    },

    saveScript: function (b,e) { 
        Ext.Msg.prompt("We have important work to do.", 
        "Enter filename.", 
        function(btn, text) {
            if (btn == 'ok'){
                yt_rpc.ExtDirectREPL.save_session({filename:text}, 
                function(f, a) {
                    if (a.result['status'] == 'SUCCESS') {
                        var alert_text = 'Saved session to ' + 
                        a.result['filename']
                        Ext.Msg.alert('Success!', alert_text);
                        var record = new logging_store.recordType(
                            {record: alert_text });
                        logging_store.add(record, number_log_records++);
                    } else {
                        Ext.Msg.alert('Always naysaying!',
                            'Failed to save to ' + 
                            a.result['filename'] + 
                            '<br>Error: ' + 
                            a.result['error']);
                    }
                });
            }
        });
    },

    downloadScript: function(b, e) {
        window.open("session.py", "_top"); 
        var record = new logging_store.recordType({
            record: 'Saved session locally.'});
        logging_store.add(record, number_log_records++);
    },

    pastebinScript: function(b, e) {
        yt_rpc.ExtDirectREPL.paste_session({}, function(f, a) {
            if (a.result['status'] == 'SUCCESS') {
                var alert_text = 'Pasted session to:<br>' + 
                a.result['site']
                var alert_text_rec = 'Pasted session to: ' + 
                a.result['site']
                Ext.Msg.alert('Pastebin', alert_text);
                var record = new logging_store.recordType(
                    {record: alert_text_rec });
                logging_store.add(record, number_log_records++);
            }
        }); 
    },

    openIRCChannel: function(b, e) {
        window.open("http://yt-project.org/irc.html", "_new");
    },

    quitReason: function(b, e) {
        task_runner.stop(heartbeat)
        yt_rpc.ExtDirectREPL.shutdown({}, function(f,a) { 
        Ext.Msg.alert("Goodbye!", "Goodbye from Reason!", function() {
        window.open("http://www.google.com/", "_top");});});
    },
});

