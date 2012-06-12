/**********************************************************************
The Parameter File widget

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

Ext.define("Reason.controller.widgets.ParameterFile", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.ParameterFileDisplay',
               'Reason.view.widgets.LevelStats',
               'Reason.store.widgets.LevelInformation'],
    templates: {
        title: "Dataset: {widget.varname}",
        createDisplay: 'widget_store.create_pf_display({varname})',
        getLevelInfo: 'widget_store["{widget.varname}"].deliver_level_stats()',
        getFieldInfo: 'widget_store["{widget.varname}"].deliver_fields()',
    },

    widgetTriggers: [
    ],

    executionTriggers: [
    ],

    viewRefs: [
        { ref:'levelStats', selector: '#levelStats'},
        { ref:'statsPanel', selector: '#statsPanel'},
    ],

    applyPayload: function(payload) {
    },

    createView: function() {
        var wd = this.payload['data'];
        examine = wd;
        this.levelDataStore = Ext.create("Reason.store.widgets.LevelInformation");
        this.levelDataStore.loadData(wd['level_stats']);
        this.levelStatsDisplay = Ext.widget("levelstats", {
            store: this.levelDataStore,
        });
        examine = this.levelStatsDisplay;
        this.dataView = Ext.widget("pfdisplay", {
             title: 'Data for ' + this.payload['varname'],
        });
        this.dataView.query("#statsPanel")[0].add(this.levelStatsDisplay);
        this.createMyRefs(this.dataView.id);
        return this.dataView;
    },

    statics: {
        widgetName: 'parameterfile',
        supportsDataObjects: false,
        supportsParameterFiles: true,
        displayName: 'Dataset Information',
        preCreation: function(obj) {
            var widget = Ext.create(this.getName());
            var cmd = widget.templateManager.applyObject(
                obj, 'createDisplay');
            reason.server.execute(cmd);
        },

    },

});
