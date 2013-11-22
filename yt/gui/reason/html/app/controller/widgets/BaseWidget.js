/**********************************************************************
Base widget class

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.require("Reason.templates.TemplateContainer");

Ext.define('Reason.controller.widgets.BaseWidget', {
    extend: 'Ext.app.Controller',
    supportsParameterFiles: false,
    supportsDataObjects: false,
    templateManager: null,
    templates: {},
    executionTriggers: [],
    widgetTriggers: [],
    keyTriggers: [],

    constructor: function() {
        this.templateManager = Ext.create(
            "Reason.templates.TemplateContainer",
            { templates: this.templates }
        );
        this.callParent(arguments);
    },

    getExecuteFunction: function(ww, templateName, isValidFn) {
        var tpl = new Ext.XTemplate(
            this.templateManager.getTemplates()[templateName]);
        var args = {};
        var control = this;
        function ev() {
            var myArgs = arguments;
            Ext.each(arguments, function(v, i) {
                args["a" + i] = myArgs[i];
            });
            args['control'] = control;
            args['widget'] = ww;
            if((isValidFn != null) && (isValidFn(arguments) == false)) {return;}
            reason.server.execute(tpl.apply(args), true);
        };
        return ev;
    },

    applyExecuteHandlers: function(ww) {
        var conf;
        function ef(id, ename, tpl, isValidFn) {
            conf = {}
            conf[ename] = this.getExecuteFunction(ww, tpl);
            ww.query(id)[0].on(conf);
        };
        Ext.each(this.executionTriggers, function(trigger) {
            ef.call(this, trigger[0], trigger[1], trigger[2], trigger[3]);
        }, this);
        Ext.each(this.widgetTriggers, function(trigger) {
            conf = {scope:this}
            conf[trigger[1]] = this[trigger[2]];
            ww.query(trigger[0])[0].on(conf);
        }, this);

        this.keyMap = new Ext.util.KeyMap({target: document});
        this.keyMap.disable();
        Ext.each(this.keyTriggers,  function(trigger) {
            trigger['fn'] = this.getExecuteFunction(ww, trigger['tpl']);
            this.keyMap.addBinding(trigger);
        }, this);
        ww.on("activate", this.activateKeyMap, this);
        ww.on("deactivate", this.deactivateKeyMap, this);
    },

    activateKeyMap: function() {
        this.keyMap.enable();
    },

    deactivateKeyMap: function() {
        this.keyMap.disable();
    },

    createMyRefs: function(varname) {
        var refs = Array(this.viewRefs.length);
        var tpl = new Ext.XTemplate("#{varname} {selector}");
        Ext.each(this.viewRefs, function(v, i, a) {
            refs[i] = {ref: v['ref'],
                       selector: tpl.apply({varname:varname,
                                            selector: v['selector']})
                      };
        });
        this.ref(refs);
        return refs;
    },

    multicast: function() {
        reason.server.multicast(this.payload['varname']);
    }

});
