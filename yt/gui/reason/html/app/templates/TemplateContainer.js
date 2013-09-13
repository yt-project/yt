/**********************************************************************
The Plot Window Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.templates.TemplateContainer", {
    constructor: function(config) {
        this.initConfig(config);
        this.callParent(arguments);
    },

    config: {
        templates : {},
    },

    applyObject: function(obj, tname) {
        var applied = {};
        var tpl;
        if (tname != null){
            tpl = new Ext.XTemplate(this.getTemplates()[tname]);
            return tpl.apply(obj);
        }
        Ext.iterate(this.getTemplates(), function(k, v, o){
            tpl = new Ext.XTemplate(v);
            applied[k] = tpl.apply(obj);
        });
        return applied;
    }
});
