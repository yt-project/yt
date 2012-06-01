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
        scrollZoom: '{widget.varname}.scroll_zoom({a1})',
        fieldChange: '{widget.varname}.set_current_field("{a1.data[\'field1\']}")',
        singleUpArrow:    '{widget.varname}.pan_rel(( 0.0, -0.1))',
        singleRightArrow: '{widget.varname}.pan_rel(( 0.1,  0.0))',
        singleDownArrow:  '{widget.varname}.pan_rel(( 0.0   0.1))',
        singleLeftArrow:  '{widget.varname}.pan_rel((-0.1,  0.0))',
        doubleUpArrow:    '{widget.varname}.pan_rel(( 0.0, -0.5))',
        doubleRightArrow: '{widget.varname}.pan_rel(( 0.5,  0.0))',
        doubleDownArrow:  '{widget.varname}.pan_rel(( 0.0   0.5))',
        doubleLeftArrow:  '{widget.varname}.pan_rel((-0.5,  0.0))',
        zoomIn10x:  '{widget.varname}.zoom(10.0)',
        zoomIn2x:   '{widget.varname}.zoom( 2.0)',
        zoomOut10x: '{widget.varname}.zoom( 0.1)',
        zoomOut2x:  '{widget.varname}.zoom( 0.5)',
        adjustTransform: '{widget.varname}.set_transform(' +
                         '{widget.varname}._current_field, ' +
                         '"{a1.data[\'field1\']}")',
        adjustColormap:  '{widget.varname}.set_cmap(' +
                         '{widget.varname}._current_field, ' +
                         '"{a1.data[\'field1\']}")',
        adjustContours:  '{widget.varname}.set_contour_info(' +
                         '"{widget.getContourField().getValue()}",' +
                         ' {widget.getNcont().getValue()},' +
                         ' {widget.getLogit().getValue():capitalize})',
        adjustVectors:   '{widget.varname}.set_vector_info(' +
                         '{widget.getVectorSkip()})',
    },

    widgetTriggers: [
        ['upload_image'],
        ['pannable_map'],
        ['clickanddrag'],
    ],

    executionTriggers: [
        ['#slider', 'changecomplete', 'scrollZoom'],
        ['#fieldSelector', 'select', 'fieldChange'],
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
        ['#transform', 'select', 'adjustTransform'],
        ['#colormap', 'select', 'adjustColormap'],
        ['#contourapply', 'click', 'adjustContours'],
        ['#vectorsapply', 'click', 'adjustVectors'],
    ],

    applyPayload: function(payload, instance) {
        this.image_panel.el.dom.src = "data:image/png;base64," + payload['image_data'];
        this.zoom_scroll.setValue(0, payload['zoom'], true);
        this.metadata_panel.update(payload['metadata_string']);
        metadata_string = payload['metadata_string'];
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
            colorbar.el.dom.src = "data:image/png;base64," +
                payload['colorbar_image'];
        }
        ticks.doLayout();
    },

    statics: {
        widgetName: 'plotwindow',
        supportsDataObjects: false,
        supportsParameterFiles: true,
        displayName: 'Do not use',
        preCreation: function(obj) {
            examine = this;
            var widget = Ext.create(this.getName())
            var ts = widget.templateManager.applyObject(obj);
            function makeProj(b, e) {
                var axis = Ext.get("axis").getValue();
                var field = Ext.get("field").getValue();
                var weight = Ext.get("weightField").getValue();
                var onmax = Ext.get("max_dens").getValue();
                reason.fireEvent("disableinput");
                yt_rpc.ExtDirectREPL.create_proj({
                        pfname: obj.varname,
                        axis: axis, field: field, weight: weight,
                        onmax: onmax},
                      Ext.emptyFn);
                Ext.WindowManager.getActive().close();
            }
            win = Ext.widget("plotwindowcreator", {title:ts.pwt, obj:obj});
            win.query("#weightField")[0].store = 
                ['None'].concat(obj.field_list);
            win.query("#field")[0].store = obj.field_list;
            /*win.show();*/
            var ww = Ext.widget("plotwindow");
            examine = ww;
            widget.applyExecuteHandlers(ww);
            Ext.ComponentQuery.query("viewport > #center-panel")[0].add(ww);
        }
    },
});
