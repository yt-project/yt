/**********************************************************************
The 3D Scene Widget

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

Ext.define("Reason.controller.widgets.Scene", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.Scene'],
    templates: {
        createScene: 'widget_store.create_scene({varname})',
        deliverGrids: 'widget_store["{widget.varname}"].deliver_gridlines()',
    },

    /* These call functions on the controller object */
    widgetTriggers: [
        ["#scenepanel", "afterrender", "setupXTK"],
    ],

    /* These call templates */
    executionTriggers: [
        ["#deliverGrids", "click", "deliverGrids"],
    ],

    /* ref: and selector: */
    viewRefs: [
        { ref:'scenePanel', selector: '#scenepanel' },
    ],

    /* key: , shift: and tpl: */
    keyTriggers: [
    ],

    applyPayload: function(payload) {
        if (payload['ptype'] == 'grid_lines') {
            this.addGridLines(payload['corners'], payload['levels']);
        } else {
            console.log("Unknown payload type received for 3D scene: " +
                        payload['ptype']);
        }
    },

    createView: function() {
        var wd = this.payload['data'];
        this.dataView = Ext.widget("scene",{
            title: wd['title'],
            varname: this.payload['varname'],
        });
        this.createMyRefs(this.dataView.id);
        this.applyExecuteHandlers(this.dataView);
        this.widgets = [];
        return this.dataView;
    },

    statics: {
        widgetName: 'scene',
        supportsDataObjects: false,
        supportsParameterFiles: true,
        displayName: '3D Plot',
        preCreation: function(obj) {
            var widget = Ext.create(this.getName());
            var cmd = widget.templateManager.applyObject(
                obj, 'createScene');
            console.log("Executing " + cmd);
            reason.server.execute(cmd);
        },
    },

    setupXTK: function() {
        var toRender = this.getScenePanel().getEl().dom.childNodes[0];
        toRender.classList.add("XTKScene");
        this.renderer = new X.renderer3D();
        this.renderer.container = toRender;
        this.renderer.init();
    },

    addGridLines: function(corners, levels) {
        var grids = new X.mesh();
        this.widgets.push(grids);
        grids.ga = "LINES";
        var p = grids.points;
        var n = grids.normals;
        var ii = 0;
        var i;
        var order1 = [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3];
        var order2 = [1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7];
        examine = corners;
        Ext.each(levels, function(level, index, allLevels) {
            for (i = 0; i < 12; i = i + 1) {
                n.add(1.0, 0.0, 0.0);
                n.add(1.0, 0.0, 0.0);
                p.add(corners[order1[i]][0][index],
                      corners[order1[i]][1][index],
                      corners[order1[i]][2][index]);
                p.add(corners[order2[i]][0][index],
                      corners[order2[i]][1][index],
                      corners[order2[i]][2][index]);
            }
        });
        this.renderer.add(grids);
        this.renderer.render();
    }

});
