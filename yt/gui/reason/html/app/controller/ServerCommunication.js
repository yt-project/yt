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
    requires: ['Reason.view.MulticastSetup', 'Reason.view.RequestsGrid'],

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
        this.requestsGrid = Ext.widget('requestsgrid');
        Ext.ComponentQuery.query("viewport > #center-panel")[0].add(
            this.requestsGrid);
    },

    execute: function(code, hide, callback) {
        var fn;
        if (callback) { fn = callback; }
        else { fn = this.returnFromRPC; }
        if (hide == null) { hide = false; }
        reason.fireEvent("disableinput");
        result_id = reason.numberOfRequests + 1;
        reason.numberOfRequests += 1;
        reason.getController("Notebook").addRequest(result_id, code);
        yt_rpc.ExtDirectREPL.execute({
            code: code,
            hide:hide,
            result_id: result_id
        }, fn);
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
            Ext.Error.raise(tpl);
        }
        reason.fireEvent("allowinput");
    },

    method: function(methodName, args, callback) {
        var m = yt_rpc.ExtDirectREPL[methodName];
        if (!m) {
            var t = "Could not identify method " + methodName;
            Ext.Msg.alert("Error", t);
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

    m: function() { return this.method(arguments); },

    multicast: function(widget_var) {
        /*
           This will be two step.  The first will be to create a window that
           allows the user to set up the MultiCast session, which will be an
           iframe loading the GAE instance.  Then we tell our local server that
           we want to multicast.
        */
        var win = Ext.widget('multicastsetup');
        setupMulticast = function() {
            var cmd = Ext.String.format(
                'widget_store.activate_multicast("{0}", "{1}", "{2}")',
                widget_var,
                win.query("#session_id")[0].getValue(),
                win.query("#session_token")[0].getValue());
            reason.server.execute(cmd);
            win.destroy();
        }
        win.query('#multicast')[0].on('click', setupMulticast);
        win.query('#cancel')[0].on('click', function(){win.destroy();});
        win.show();
    },

});

