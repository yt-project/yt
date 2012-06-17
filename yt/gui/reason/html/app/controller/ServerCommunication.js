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

Ext.define('Reason.controller.ServerCommunication', {
    extend: 'Ext.app.Controller',
    stores: ['Payloads'],
    views: ['PayloadGrid'],

    init: function() {
        this.application.addListener({
            payloadreceived: {fn: this.handlePayload, scope: this},
            stopheartbeat: {fn: this.stopHeartbeat, scope: this},
            enabledebug: {fn: this.enableDebug, scope: this},
        })
        /* We also use this as our heartbeat */
        this.taskRunner = new Ext.util.TaskRunner();
        heartbeatRequest = false;
        this.heartbeat = this.taskRunner.start(
            {run: this.heartbeatCall,
             interval: 250});
        this.dataObjectBeat = this.taskRunner.start(
            {run: this.dataObjectsCall,
             interval: 5000});
        this.callParent(arguments);
    },

    handlePayload: function(payload) {
        if ((this.payloadGrid) && (payload['type'] != 'logentry')) {
            var wv = payload['varname'];
            if (payload['type'] == 'widget_payload') {
                wv = payload['widget_id'];
            }
            this.getPayloadsStore().add({
                payloadtype: payload['type'],
                widgettype: payload['widget_type'],
                varname: wv,
            });
        }
        this.application.fireEvent('payload' + payload['type'], payload);
    },

    dataObjectsCall: function() {
        yt_rpc.ExtDirectParameterFileList.get_list_of_pfs({}, 
            function(f, a) {
                if (f == null) { return; }
                reason.fireEvent("newdataobjects", f);
        });
    },

    heartbeatCall: function() {
        if (heartbeatRequest == true) return;
        heartbeatRequest = true;
        console.log("Sending heartbeat");
        yt_rpc.ExtDirectREPL.heartbeat(
            {}, function(f, a) {
                heartbeatRequest = false;
                if (f != null) { 
                    Ext.each(f, function(payload, index) {
                            reason.fireEvent("payloadreceived", payload);
                    });
                }
                /*else if (!a.status) {
                    reason.fireEvent("stopheartbeat");
                    Ext.Msg.alert("Error", "Error talking to server.  Shutting down.");
                }*/
                return true;
            }
        );
    },

    stopHeartbeat: function() {
        this.taskRunner.stopAll();
    },

    enableDebug: function() {
        this.payloadGrid = Ext.widget("payloadgrid");
        Ext.ComponentQuery.query("viewport > #center-panel")[0].add(
            this.payloadGrid);
    },

    execute: function(code, hide, callback) {
        var fn;
        if (callback) { fn = callback; }
        else { fn = this.returnFromRPC; }
        if (hide == null) { hide = false; }
        reason.fireEvent("disableinput");
        yt_rpc.ExtDirectREPL.execute({code: code, hide:hide}, fn);
    },

    returnFromRPC: function(result, e) {
        if(!e.status) {
            var tpl = new Ext.XTemplate(
                'RPC Error: {message}; {action}, {method}');
            var trans = e.getTransaction();
            tpl = tpl.apply({message: e.message,
                             action: trans.action,
                             method: trans.method});
            Ext.Msg.alert("Error", tpl);
            examine = {result: result, e: e};
            Ext.Error.raise(tpl);
        }
        reason.fireEvent("allowinput");
    },

    method: function(methodName, args, callback) {
        var m = yt_rpc.ExtDirectREPL[methodName];
        if (!m) {
            var t = "Could not identify method " + methodName;
            Ext.Msg.alert("Error", t);
            examine = {result: result, e: e};
            Ext.Error.raise(t);
        }
        var fn;
        if (callback) { fn = callback; }
        else {
            this.application.fireEvent("disableinput");
            fn = this.returnFromRPC;
        }
        m(args, fn);
    },

    m: function() { return this.method(arguments); }

});

