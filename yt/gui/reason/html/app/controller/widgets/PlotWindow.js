/**********************************************************************
The Plot Window Widget

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

Ext.define("Reason.controller.widgets.PlotWindow", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.PlotWindowCreator',
               'Reason.view.widgets.PlotWindow'],
    templates: {
        pwt: 'Projection Details for {name}',
        swt: 'Slice Details for {name}',
        refresh: 'widget_store["{widget.varname}"].refresh()',
        createSlice: 'widget_store.create_slice({pfname}, [{center}], "{axis}",' + 
                     '"{field}", {onmax:capitalize})',
        createProj: 'widget_store.create_proj({pfname}, "{axis}",' + 
                     '"{field}", "{weight}")',
        scrollZoom: 'widget_store["{widget.varname}"].scroll_zoom({a1})',
        fieldChange: 'widget_store["{widget.varname}"].set_current_field("{a1}")',
        singleUpArrow:    'widget_store["{widget.varname}"].pan_rel(( 0.0, -0.1))',
        singleRightArrow: 'widget_store["{widget.varname}"].pan_rel(( 0.1,  0.0))',
        singleDownArrow:  'widget_store["{widget.varname}"].pan_rel(( 0.0,  0.1))',
        singleLeftArrow:  'widget_store["{widget.varname}"].pan_rel((-0.1,  0.0))',
        doubleUpArrow:    'widget_store["{widget.varname}"].pan_rel(( 0.0, -0.5))',
        doubleRightArrow: 'widget_store["{widget.varname}"].pan_rel(( 0.5,  0.0))',
        doubleDownArrow:  'widget_store["{widget.varname}"].pan_rel(( 0.0,  0.5))',
        doubleLeftArrow:  'widget_store["{widget.varname}"].pan_rel((-0.5,  0.0))',
        zoomIn10x:  'widget_store["{widget.varname}"].zoom(10.0)',
        zoomIn2x:   'widget_store["{widget.varname}"].zoom( 2.0)',
        zoomOut10x: 'widget_store["{widget.varname}"].zoom( 0.1)',
        zoomOut2x:  'widget_store["{widget.varname}"].zoom( 0.5)',
        adjustTransform: 'widget_store["{widget.varname}"].set_transform(' +
                         'widget_store["{widget.varname}"]._current_field, ' +
                         '"{a1}")',
        adjustColormap:  'widget_store["{widget.varname}"].set_cmap(' +
                         'widget_store["{widget.varname}"]._current_field, ' +
                         '"{a1}")',
        adjustContours:  'widget_store["{widget.varname}"].set_contour_info(' +
                         '"{[values.control.getContourField().getValue()]}",' +
                         ' {[values.control.getNumContours().getValue()]},' +
                         ' {[(Ext.util.Format.capitalize(""+values.control.getLogContours().getValue()))]})',
        adjustVectors:   'widget_store["{widget.varname}"].set_vector_info(' +
                         '{[values.control.getVectorSkip().getValue()]})',
        recenterImage:   'widget_store["{widget.varname}"].image_recenter(' +
                         '{x}, {y}, {w}, {h})',
        dragImage:       'widget_store["{widget.varname}"].pan_rel(({rel_x}, {rel_y}))',
        createPannableMap: 'widget_store.create_mapview("{widget.varname}")',
    },

    widgetTriggers: [
        ['#uploadimage', 'click', 'uploadImage'],
        ['#imagepanel', 'afterrender', 'setupClickImage'],
        ['#multicast', 'click', 'multicast'],
    ],

    executionTriggers: [
        ['#zoomSlider', 'changecomplete', 'scrollZoom'],
        ['#fieldSelector', 'change', 'fieldChange'],
        ['#singleuparrow', 'click', 'singleUpArrow'],
        ['#singlerightarrow', 'click', 'singleRightArrow'],
        ['#singledownarrow', 'click', 'singleDownArrow'],
        ['#singleleftarrow', 'click', 'singleLeftArrow'],
        ['#doubleuparrow', 'click', 'doubleUpArrow'],
        ['#doublerightarrow', 'click', 'doubleRightArrow'],
        ['#doubledownarrow', 'click', 'doubleDownArrow'],
        ['#doubleleftarrow', 'click', 'doubleLeftArrow'],
        ['#zoomin10x', 'click', 'zoomIn10x'],
        ['#zoomin2x', 'click', 'zoomIn2x'],
        ['#zoomout10x', 'click', 'zoomOut10x'],
        ['#zoomout2x', 'click', 'zoomOut2x'],
        ['#transform', 'change', 'adjustTransform'],
        ['#colormap', 'change', 'adjustColormap'],
        ['#contourapply', 'click', 'adjustContours'],
        ['#vectorsapply', 'click', 'adjustVectors'],
        ['#pannablemap', 'click', 'createPannableMap'],
    ],

    viewRefs: [
        { ref: 'colorbar', selector: '#colorbar'},
        { ref: 'image', selector: '#imagepanel'},
        { ref: 'fieldSelector', selector: '#fieldSelector'},
        { ref: 'transform', selector: '#transform'},
        { ref: 'contourField', selector: '#contourfield'},
        { ref: 'numContours', selector: '#ncont'},
        { ref: 'logContours', selector: '#logit'},
        { ref: 'zoomSlider', selector: '#zoomSlider'},
        { ref: 'metadataString', selector: '#metadataString'},
        { ref: 'ticks', selector: '#ticks'},
        { ref: 'vectorSkip', selector: '#skip'},
    ],

    keyTriggers: [
        { key: 'z', shift: false, tpl: "zoomIn2x" },
        { key: 'Z', shift: true, tpl: "zoomIn10x" },
        { key: 'x', shift: false, tpl: "zoomOut2x" },
        { key: 'X', shift: true, tpl: "zoomOut10x" },
        { key: 'k', shift: false, tpl: "singleUpArrow" },
        { key: 'j', shift: false, tpl: "singleDownArrow" },
        { key: 'h', shift: false, tpl: "singleLeftArrow" },
        { key: 'l', shift: false, tpl: "singleRightArrow" },
        { key: 'K', shift: true, tpl: "doubleUpArrow" },
        { key: 'J', shift: true, tpl: "doubleDownArrow" },
        { key: 'H', shift: true, tpl: "doubleLeftArrow" },
        { key: 'L', shift: true, tpl: "doubleRightArrow" },
    ],

    applyPayload: function(payload) {
        this.getImage().getEl().dom.src = 
            "data:image/png;base64," + payload['image_data'];
        this.getZoomSlider().setValue(0, payload['zoom'], true);
        this.getMetadataString().update(payload['metadata_string']);
        this.getMetadataString().mds = payload['metadata_string'];
        var ticks = this.getTicks();
        ticks.removeAll();
        Ext.each(payload['ticks'], function(tick, index) {
            ticks.add({xtype:'panel',
                       width: 10, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:0, y: 10 + tick[0]});
            ticks.add({xtype:'panel',
                       width: 90, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px;',
                       html: '' + tick[2] + '',
                       x:12, y: 4 + tick[0]});
        });
        if (payload['colorbar_image'] != null) {
            this.getColorbar().getEl("img").dom.src =
                "data:image/png;base64," + payload['colorbar_image'];
        }
        ticks.doLayout();
    },

    createView: function() {
        var wd = this.payload['data'];
        this.plotWindowView = Ext.widget("plotwindow",{
            varname : this.payload['varname'],
            title: wd['title'],
        });
        var newRefs = this.createMyRefs(this.plotWindowView.id);
        var refresh = this.templateManager.applyObject(
            {widget: this.plotWindowView}, 'refresh');
        reason.server.execute(refresh, false);
        this.getColorbar().src = "data:image/png;base64," + wd['colorbar'];
        this.fieldStore = Ext.create("Reason.store.Fields")
        this.fieldStore.loadData(wd['fields']);
        this.getFieldSelector().bindStore(this.fieldStore);
        this.getFieldSelector().setValue(wd['initial_field']);
        this.getTransform().setValue(wd['initial_transform']);
        this.getContourField().bindStore(this.fieldStore);
        this.getContourField().setValue(wd['initial_field']);
        this.applyExecuteHandlers(this.plotWindowView);
        return this.plotWindowView;
    },

    createMulticastView: function() {
        this.plotWindowView = Ext.widget("plotwindow", {
            varname: null,
            title: "Multicast View",
        });
        Ext.each(this.widgetTriggers, function(trigger, index, all) {
            this.plotWindowView.query(trigger[0])[0].disable();
        });
        Ext.each(this.executionTriggers, function(trigger, index, all) {
            this.plotWindowView.query(trigger[0])[0].disable();
        });
    },

    statics: {
        widgetName: 'plot_window',
        supportsDataObjects: false,
        supportsParameterFiles: true,
        displayName: '2D Plot',
        preCreation: function(obj) {
            var widget = Ext.create(this.getName());
            function makePlot(b, e) {
                var conf = {
                    pfname: obj.varname,
                    axis: win.query("#axis")[0].getValue(),
                    field: win.query("#field")[0].getValue(),
                    onmax: "" + win.query("#maxDens")[0].getValue(),
                };
                var method = 'createSlice';
                if (win.query("#plotType")[0].getValue() == 'Projection') {
                    method = 'createProj';
                    conf['weight'] = win.query("#weightField")[0].getValue();
                } else {
                  conf['center'] = [win.query("#slice_x_center")[0].getValue(),
                                    win.query("#slice_y_center")[0].getValue(),
                                    win.query("#slice_z_center")[0].getValue()];
                }
                var cmd = widget.templateManager.applyObject(conf, method);
                reason.server.execute(cmd);
                win.destroy();
            }
            function togglePlotType(b, e) {
                var plotType = win.query("#plotType")[0].getValue();
                if (plotType == 'Projection') {
                    win.query("#weightField")[0].enable();
                    win.query("#maxDens")[0].disable();
                    win.query("#slice_x_center")[0].disable();
                    win.query("#slice_y_center")[0].disable();
                    win.query("#slice_z_center")[0].disable();
                } else {
                    win.query("#weightField")[0].disable();
                    win.query("#maxDens")[0].enable();
                    win.query("#slice_x_center")[0].enable();
                    win.query("#slice_y_center")[0].enable();
                    win.query("#slice_z_center")[0].enable();
                }
            }
            function toggleMaxDens(checkbox, checked) {
                var plotType = win.query("#plotType")[0].getValue();
                if (plotType == "Projection") { return; }
                if (checked == true) {
                    win.query("#slice_x_center")[0].disable();
                    win.query("#slice_y_center")[0].disable();
                    win.query("#slice_z_center")[0].disable();
                } else {
                    win.query("#slice_x_center")[0].enable();
                    win.query("#slice_y_center")[0].enable();
                    win.query("#slice_z_center")[0].enable();
                }
            }
            var title = widget.templateManager.applyObject(obj, 'pwt');
            win = Ext.widget("plotwindowcreator", {title:title, obj:obj});
            console.log(obj.field_list);
            var field_list = [];
            Ext.each(obj.field_list, function(f, i, af) {
                field_list.push(f.text);
            });
            win.query("#weightField")[0].store = ['None'].concat(field_list);
            win.query("#field")[0].store = field_list;
            win.query("#create")[0].on('click', makePlot);
            win.query("#cancel")[0].on('click', function(){win.destroy();});
            win.query("#maxDens")[0].on('change', toggleMaxDens);
            win.query("#plotType")[0].on('change', togglePlotType);
            togglePlotType();
            toggleMaxDens();
            win.show();
            /* Note that in this case, our instance of 'widget', which is this
               class, is not long-lived.  It dies after the window is
               destroyed. */
        },

    },

    uploadImage: function() {
        var imageData = this.getImage().getEl().dom.src;
        var mds = this.getMetadataString().mds;
        yt_rpc.ExtDirectREPL.upload_image(
            {image_data:imageData, caption:mds},
            function(rv) {
                var alert_text;
                if(rv['uploaded'] == false) {
                    alert_text = "Failure uploading image!";
                } else {
                    alert_text = "Uploaded to " +
                            rv['upload']['links']['imgur_page'];
                }
                reason.fireEvent("logentry", alert_text);
                Ext.Msg.alert('imgur.com', alert_text);
            });
    },

    setupClickImage: function(c) {
        var controller = this;
        var dragStartPos = {x:-1, y:-1};
        c.el.on('click',  function(e) {
            if (e.ctrlKey == false) return;
            var xy = e.getXY();
            var args = {widget: controller.plotWindowView,
                        x: xy[0] - controller.getImage().getEl().dom.x,
                        y: xy[1] - controller.getImage().getEl().dom.y,
                        w: controller.getImage().width,
                        h: controller.getImage().height};
            var code = controller.templateManager.applyObject(
                        args, "recenterImage");
            reason.server.execute(code, true);
        });
    },
});
