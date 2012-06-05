/**********************************************************************
The Plot Window Widget

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

// shim layer with setTimeout fallback
Ext.define("Reason.controller.widgets.GridDataViewer", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.GridDataViewer'],

    templates: {
        title: "Grid Data for {widget.varname}",
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

        examine = {wd:wd, tt:this};

        this.createMyRefs(this.gridDataView.id);
        return this.gridDataView;
    },

    statics: {
        widgetName: "grid_data",
        displayName: "Grid Data Viewer",
        supportsDataObjects: false,
        supportsParameterFiles: true,
        preCreation: function(obj) {
            reason.server.method("create_grid_dataview", {pfname:obj.varname});
        },
    },
});
