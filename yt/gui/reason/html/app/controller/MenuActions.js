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
    views: ['MainMenu'],

    init: function() {
        this.control({
            '#openfile': { click: this.openFile },
            '#savescript': { click: this.saveScript},
            '#downloadscript': { click: this.downloadScript },
            '#pastebinscript': { click: this.pastebinScript },
            '#ytchat': {click: this.openIRCChannel },
            '#quit': {click: this.quitReason },
            '#enabledebug': {click: this.enableDebug },
        });
        this.callParent(arguments);
    },

    saveScript: function (b,e) { 
        /* This function gets called on success */
        var controller = this;
        function handleResponse(f, a) {
            if (a.result['status'] == 'SUCCESS') {
                var alert_text = 'Saved session to ' + a.result['filename'];
                Ext.Msg.alert('Success!', alert_text);
                controller.application.fireEvent("logentry", alert_text);
            } else {
                Ext.Msg.alert('Always naysaying!',
                  'Failed to save to ' + a.result['filename'] + '<br>Error: ' +
                  a.result['error']);
           }
        };
        /* Now we prompt */
        Ext.Msg.prompt("We have important work to do.", "Enter filename.", 
            function(btn, text) {
                if (btn == 'ok') {
                    reason.server.method('save_session',
                        {filename:text}, handleResponse);
                }
            }
        );
    },

    openFile: function(b, e) {
        alert("Not yet implemented.");
    },

    downloadScript: function(b, e) {
        window.open("session.py", "_top"); 
        this.application.fireEvent("logentry", 'Saved session locally.')
    },

    pastebinScript: function(b, e) {
        reason.server('paste_session', {}, function(f, a) {
            if (a.result['status'] == 'SUCCESS') {
                var alert_text = 'Pasted session to:<br>' + 
                a.result['site']
                var alert_text_rec = 'Pasted session to: ' + 
                a.result['site']
                Ext.Msg.alert('Pastebin', alert_text);
                this.application.fireEvent("logentry", alert_text_rec);
            }
        }); 
    },

    openIRCChannel: function(b, e) {
        window.open("http://yt-project.org/irc.html", "_new");
    },

    quitReason: function(b, e) {
        this.application.fireEvent("stopheartbeat");
        reason.server('shutdown', {}, function(f, a) {
        Ext.Msg.alert("Goodbye!", "Goodbye from Reason!", function() {
        window.open("http://www.google.com/", "_top");});});
    },

    enableDebug: function() {
        this.application.fireEvent("enabledebug");
    },
});

