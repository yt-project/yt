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
               'Reason.view.widgets.FieldPanel',
               'Reason.store.widgets.LevelInformation'],
    templates: {
        title: "Dataset: {widget.varname}",
        createDisplay: 'widget_store.create_pf_display({varname})',
        retrieveField: 'widget_store["{widget.varname}"].deliver_field("{a1}")',
        fieldSource: '<pre>{field_source}</pre>',
        pfParams: '<table class="pftable">' +
                  '<tr><th>Parameter</th><th>Value</th></tr>' +
                  '<tr><td>Output Hash</td>' + 
                  '<td>{output_hash}</td></tr>' +
                  '<tr><td>Dimensionality</td>' +
                  '<td>{dimensionality}</td></tr>' +
                  '<tr><td>Refine by</td>' +
                  '<td>{refine_by}</td></tr>' +
                  '<tr><td>Domain Dimensions</td>' +
                  '<td>{domain_dimensions}</td></tr>' +
                  '<tr><td>Cosmological Simulation</td>' +
                  '<td>{cosmological_simulation}</td></tr>' +
                  '<tr><td>Current Redshift</td>' +
                  '<td>{current_redshift}</td></tr>' +
                  '<tr><td>Omega Matter</td>' +
                  '<td>{omega_matter}</td></tr>' +
                  '<tr><td>Omega Lambda</td>' +
                  '<td>{omega_lambda}</td></tr>' +
                  '<tr><td>Hubble Constant</td>' +
                  '<td>{hubble_constant}</td></tr>' +
                  '<tr><td>Current (sim) Time</td>' +
                  '<td>{current_time}</td></tr>' +
                  '<tr><td>Domain Left Edge</td>' +
                  '<td>{domain_left_edge}</td></tr>' +
                  '<tr><td>Domain Right Edge</td>' +
                  '<td>{domain_right_edge}</td></tr>' +
                  '</table>',
    },

    widgetTriggers: [
    ],

    executionTriggers: [
        ['#fieldSelector', 'change', 'retrieveField'],
    ],

    viewRefs: [
        { ref:'levelStats', selector: '#levelStats'},
        { ref:'statsPanel', selector: '#statsPanel'},
        { ref:'parametersPanel', selector: '#pfParams'},
        { ref:'fieldPanel', selector: '#fieldPanel'},
        { ref:'fieldSourcePanel', selector: '#fieldSourcePanel'},
        { ref:'widgetPanel', selector: '#widgetpanel'},
    ],

    applyPayload: function(payload) {
        if (payload['ptype'] == 'field_info') {
            var source = this.templateManager.applyObject(
                payload, 'fieldSource')
            this.getFieldSourcePanel().update(source);
        }
    },

    createView: function() {
        var wd = this.payload['data'];
        examine = wd;
        this.levelDataStore = Ext.create("Reason.store.widgets.LevelInformation");
        this.levelDataStore.loadData(wd['level_stats']);
        this.levelStatsDisplay = Ext.widget("levelstats", {
            store: this.levelDataStore,
        });
        this.dataView = Ext.widget("pfdisplay", {
             title: 'Data for ' + this.payload['varname'],
             varname: this.payload['varname'],
        });
        this.fieldDisplay = Ext.widget("fieldpanel");
        fieldSelector = this.fieldDisplay.query("#fieldSelector")[0]
        this.fieldStore = Ext.create("Reason.store.Fields")
        this.fieldStore.loadData(wd['fields']);
        fieldSelector.bindStore(this.fieldStore);

        this.dataView.query("#fieldPanel")[0].add(this.fieldDisplay);
        this.dataView.query("#statsPanel")[0].add(this.levelStatsDisplay);
        this.createMyRefs(this.dataView.id);

        var pfString = this.templateManager.applyObject(
                wd['pf_info'], 'pfParams');
        this.getParametersPanel().update(pfString);
        this.applyExecuteHandlers(this.dataView);
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
