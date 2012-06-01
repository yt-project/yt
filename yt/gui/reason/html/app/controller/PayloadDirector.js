/**********************************************************************
The Payload handling facility for Reason

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

var heartbeatRequest = false;

Ext.define('Reason.controller.PayloadDirector', {
    extend: 'Ext.app.Controller',

    init: function() {
        this.application.addListener({
            payloadreceived: {fn: this.handlePayload, scope: this},
            stopheartbeat: {fn: this.stopHeartbeat, scope: this},
        })
        /* We also use this as our heartbeat */
        this.taskRunner = new Ext.util.TaskRunner();
        heartbeatRequest = false;
        this.heartbeat = this.taskRunner.start(
            {run: this.heartbeatCall,
             interval: 250});
        this.heartbeat = this.taskRunner.start(
            {run: this.dataObjectsCall,
             interval: 5000});
        this.callParent(arguments);
    },
    handlePayload: function(payload) {
        this.application.fireEvent('payload' + payload['type'], payload);
    },

    dataObjectsCall: function() {
        yt_rpc.ExtDirectParameterFileList.get_list_of_pfs({}, 
            function(f, a) {
                if (f == null) { alert("Error!"); return; }
                reason.fireEvent("newdataobjects", f);
        });
    },

    heartbeatCall: function() {
        if (heartbeatRequest == true) return;
        heartbeatRequest = true;
        yt_rpc.ExtDirectREPL.heartbeat(
            {}, function(f, a) {
                heartbeatRequest = false;
                if (f != null) { 
                    Ext.each(f, function(payload, index) {
                            reason.fireEvent("payloadreceived", payload);
                    });
                }
                return true;
            }
        );
    },

    stopHeartbeat: function() {
        this.taskRunner.stopAll();
    },
});

