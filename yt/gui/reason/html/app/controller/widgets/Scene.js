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
    requires: ['Reason.view.widgets.Scene',
               'Reason.view.widgets.IsocontourCreator',
               'Reason.store.widgets.CameraKeyFrames',
               'Reason.store.widgets.CameraPathElements',
               'Reason.store.widgets.SceneWidgets',
               'Ext.ux.CheckColumn'],
    templates: {
        createScene: 'widget_store.create_scene({varname})',
        deliverGrids: 'widget_store["{widget.varname}"].deliver_gridlines()',
        createIsocontour: 'widget_store["{varname}"].deliver_isocontour(' +
                          '"{field}", {value}, {relValue:capitalize})',
    },

    /* These call functions on the controller object */
    widgetTriggers: [
        ["#scenepanel", "afterrender", "setupXTK"],
        ["#addKeyframe", "click", "addKeyframe"],
        ["#renderPath", "click", "renderPath"],
        ["#keyframeview", "select", "shiftToKeyframe"],
        ["#widgetEnabled", "checkchange", "toggleWidgetEnabled"],
        ["#addIsocontour", "click", "createIsocontour"],
        ["#cameraPathSlider", "change", "updateCameraPosition"],
    ],

    /* These call templates */
    executionTriggers: [
    ],

    /* ref: and selector: */
    viewRefs: [
        { ref:'scenePanel', selector: '#scenepanel' },
        { ref:'keyFrameView', selector: '#keyframeview' },
        { ref:'widgetPanel', selector: '#widgetlist'},
        { ref:'cameraPathSlider', selector: '#cameraPathSlider'},
    ],

    /* key: , shift: and tpl: */
    keyTriggers: [
    ],

    applyPayload: function(payload) {
        if (payload['ptype'] == 'grid_lines') {
            payload['corners'] = new Float64Array(payload['corners']);
            payload['levels'] = new Int32Array(payload['levels']);
            this.addGridLines(payload['corners'], payload['levels'], payload['max_level']);
        } else if (payload['ptype'] == 'isocontour') {
            payload['vert'] = new Float64Array(payload['vert']);
            payload['normals'] = new Float64Array(payload['normals']);
            this.addIsocontour(payload['vert'], payload['normals']);
        } else if (payload['ptype'] == 'camerapath') {
            this.updateCameraPathElements(payload['data']);
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
        this.keyFrames = Ext.create("Reason.store.widgets.CameraKeyFrames");
        this.getKeyFrameView().bindStore(this.keyFrames);
        this.pathElements = Ext.create("Reason.store.widgets.CameraPathElements");
        this.widgets = Ext.create("Reason.store.widgets.SceneWidgets");
        this.getWidgetPanel().bindStore(this.widgets);
        this.fieldStore = Ext.create("Reason.store.Fields")
        this.fieldStore.loadData(wd['fields']);
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
        this.renderer.interactor.config.KEYBOARD_ENABLED = false;
        this.renderer.interactor.init();
        var cmd = this.templateManager.applyObject(
            {widget: this.dataView}, 'deliverGrids');
        reason.server.execute(cmd);
    },

    addGridLines: function(corners, levels, maxLevel) {
        var i, g, n, p, offset, ind;
        var order1 = [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3];
        var order2 = [1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7];
        var nv = levels.length;
        gw = [];
        for (i = 0; i < maxLevel + 1; i = i + 1) {
            var grids = new X.mesh();
            gw.push(grids);
            this.widgets.add({
              name: "Grid Level " + i,
              type: 'grids',
              widget: grids,
              enabled: true,
            });
            grids.ga = "LINES";
        }
        examine = {n: n, p: p, corners: corners};
        var i0, i1, i2;
        Ext.each(levels, function(level, index, allLevels) {
            p = gw[level].points;
            n = gw[level].normals;
            for (i = 0; i < 12; i = i + 1) {
                n.add(1.0, 0.0, 0.0);
                n.add(1.0, 0.0, 0.0);
                p.add(corners[(((order1[i] * 3 + 0)*nv)+index)],
                      corners[(((order1[i] * 3 + 1)*nv)+index)],
                      corners[(((order1[i] * 3 + 2)*nv)+index)]);
                p.add(corners[(((order2[i] * 3 + 0)*nv)+index)],
                      corners[(((order2[i] * 3 + 1)*nv)+index)],
                      corners[(((order2[i] * 3 + 2)*nv)+index)]);
            }
        });
        for (i = 0; i < maxLevel + 1; i = i + 1) {
            this.renderer.add(gw[i]);
        }
        this.renderer.render();
    },

    createIsocontour: function() {
        var win; 
        var controller = this;
        /* field , value */
        function callExtract(b, e) {
            var conf = {
                varname: controller.dataView.varname,
                field: win.query("#field")[0].getValue(),
                value: win.query("#value")[0].getValue(),
                relValue: "" + win.query("#relValue")[0].getValue(),
            };
            cmd = controller.templateManager.applyObject(
                    conf, "createIsocontour");
            reason.server.execute(cmd);
            win.destroy();
        }
        win = Ext.widget("isocontourcreator");
        win.query("#field")[0].bindStore(this.fieldStore);
        win.query("#extract")[0].on('click', callExtract);
        win.query("#cancel")[0].on('click', function(){win.destroy();});
        win.show();
    },

    addIsocontour: function(vertices, normals) {
        console.log("Adding isocontours ...");
        var i, g, n, p;
        var nv = vertices.length/3;
        var last = 0;
        var surf = new X.mesh();
        this.widgets.add({
          name: "Isocontour",
          type: 'isocontour',
          widget: surf,
          enabled: true,
        });
        surf.ga = "TRIANGLES";
        p = surf.points;
        n = surf.normals;

        for (index = 0; index < nv; index = index + 1) {
            p.add(vertices[index * 3 + 0],
                  vertices[index * 3 + 1],
                  vertices[index * 3 + 2]);
            n.add(normals[index * 3 + 0],
                  normals[index * 3 + 1],
                  normals[index * 3 + 2]);
        }
        surf.color = [1.0, 0.0, 0.0];
        this.renderer.add(surf);
        this.renderer.render();
    },

    addKeyframe: function() {
        this.getCameraPathSlider().setValue(0);
        this.getCameraPathSlider().disable();
        var v = this.renderer.camera.view;
        var va = v.toArray();
        this.keyFrames.add({
            time: 1.0,
            view: v,
            pos_x: -va[0][3],
            pos_y: -va[1][3],
            pos_z: -va[2][3],
        });
    },

    shiftToKeyframe: function(rowModel, record, index) {
        this.renderer.camera.view = record.data.view;
        this.renderer.render();
    },

    toggleWidgetEnabled: function(column, rowIndex, enabled) {
        /* We have to get the store, then the widget ... */
        var rec = this.widgets.data.items[rowIndex];
        rec.data.widget.visible = enabled;
        this.renderer.render();
    },

    renderPath: function() {
        var t = new Ext.XTemplate("[[{0}, {1}, {2}, {3}], ",
                                  " [{4}, {5}, {6}, {7}], ",
                                  " [{8}, {9}, {10}, {11}], ",
                                  " [{12}, {13}, {14}, {15}]],\n");
        var cmdt = new Ext.XTemplate("widget_store['{0}'].render_path(\n",
                                     "[{1}]\n,[{2}], 101)");
        var path = "";
        var times = "";
        Ext.each(this.keyFrames.data.items, function(rec, ind, all) {
            times = times + rec.data.time + ", ";
            path = path + t.apply(rec.data.view.flatten());
        });
        var cmd = cmdt.apply([this.dataView.varname, path, times]);
        reason.server.execute(cmd);
    },

    updateCameraPathElements: function(elements) {
        var cpe = this.pathElements;
        var i;
        cpe.removeAll();
        for (i = 0; i < elements[0].length; i = i + 1) {
            cpe.add({position: elements[0][i],
                     focus: elements[1][i],
                     up: elements[2][i]});
        }
        v = this.getCameraPathSlider().enable();
    },

    updateCameraPosition: function(b, e) {
        v = this.getCameraPathSlider().getValue();
        console.log(v);
        rec = this.pathElements.data.items[v].data;
        this.renderer.camera.position = rec.position;
        this.renderer.camera.focus = rec.focus;
        this.renderer.camera.up = rec.up;
        this.renderer.render();
    },

});
