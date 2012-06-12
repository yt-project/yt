/**********************************************************************
Base widget class

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

Ext.require("Reason.templates.TemplateContainer");

Ext.define('Reason.controller.widgets.BaseWidget', {
    extend: 'Ext.app.Controller',
    supportsParameterFiles: false,
    supportsDataObjects: false,
    templateManager: null,
    templates: {},
    executionTriggers: [],
    widgetTriggers: [],

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
            console.log("Inside ... " + templateName);
            var myArgs = arguments;
            Ext.each(arguments, function(v, i) {
                args["a" + i] = myArgs[i];
            });
            args['control'] = control;
            args['widget'] = ww;
            examine = {args:args, arg2: arguments};
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
            /*console.log(trigger[0] + " " + trigger[1] + " " + trigger[2], trigger[3]);*/
            ef.call(this, trigger[0], trigger[1], trigger[2], trigger[3]);
        }, this);
        Ext.each(this.widgetTriggers, function(trigger) {
            console.log(trigger[0], trigger[1], trigger[2]);
            conf = {scope:this}
            conf[trigger[1]] = this[trigger[2]];
            ww.query(trigger[0])[0].on(conf);
        }, this);
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

});
