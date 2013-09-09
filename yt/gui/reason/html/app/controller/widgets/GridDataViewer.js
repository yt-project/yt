/**********************************************************************
The Plot Window Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

// shim layer with setTimeout fallback
Ext.define("Reason.controller.widgets.GridDataViewer", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.GridDataViewer'],

    templates: {
        title: "Grid Data for {widget.varname}",
        createGridViewer: 'widget_store.create_grid_dataview({varname})',
    },

    widgetTriggers: [

    ],

    executionTriggers: [

    ],

    viewRefs: [
    ],

    applyPayload: function(payload) {
        return;
    },

    createView: function() {
        var wd = this.payload['data'];
        this.dataStore = Ext.create("Reason.store.widgets.GridData");
        this.dataStore.loadData(wd['gridvals']);
        this.gridDataView = Ext.widget("griddataviewer", {
             store: this.dataStore,
             title: 'Grid Data for ' + this.payload['varname'],
        });

        this.createMyRefs(this.gridDataView.id);
        return this.gridDataView;
    },

    statics: {
        widgetName: "grid_data",
        displayName: "Grid Data Viewer",
        supportsDataObjects: false,
        supportsParameterFiles: true,
        preCreation: function(obj) {
            var widget = Ext.create(this.getName());
            var cmd = widget.templateManager.applyObject(
                obj, 'createGridViewer');
            reason.server.execute(cmd);
        },
    },
});
