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


Ext.define("Reason.controller.widgets.PhasePlot", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.PhasePlot'],

    viewRefs: [
        { ref: 'yTicks', selector: '#y_ticks'},
        { ref: 'xTicks', selector: '#x_ticks'},
        { ref: 'colorTicks', selector: 'cb_ticks'},
        { ref: 'colorbar', selector: '#colorbar'},
        { ref: 'image', selector: '#imagepanel'},
        { ref: 'metadataString', selector: '#metadataString'},
    ],
    
    applyPayload: function(payload) {
        this.getImage().getEl().dom.src = 
            "data:image/png;base64," + payload['image_data'];
        this.getMetadataString().update(payload['metadata_string']);
        this.getMetaDataString().mds = payload['metadata_string'];
        ticks.removeAll();
        this.getColorbar.getEl().dom.src=
            "data:image/png;base64," + payload['cbar']['cmap_image'];
        var YTicks = this.getYTicks();
        Ext.each(payload['yax']['ticks'], function(tick, index) {
            YTicks.add({xtype:'panel',
                       width: 20, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: right;' +
                              'padding-right: 5px;',
                       html: ' ' + tick[2] + ' ',
                       x:0, y: tick[0]-6});
            YTicks.add({xtype:'panel',
                       width: 20, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:20, y: tick[0]});
        });
        YTicks.doLayout();
        var XTicks = this.getXTicks();
        Ext.each(payload['xax']['ticks'], function(tick, index) {
            XTicks.add({xtype:'panel',
                       width: 1, height:20,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:(400 - tick[0]) + 10, y: 0});
            XTicks.add({xtype:'panel',
                       width: 20, height:20,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: center;',
                       html: ' ' + tick[2] + ' ',
                       x: (400 - tick[0]), y: 20});
        });
        XTicks.doLayout();
        var colorTicks = this.getColorTicks();
        Ext.each(payload['cbar']['ticks'], function(tick, index) {
            colorTicks.add({xtype:'panel',
                       width: 10, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:0, y: tick[0]});
            colorTicks.add({xtype:'panel',
                       width: 30, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px;',
                       html: ' ' + tick[2] + ' ',
                       x:18, y: tick[0]-6});
        });
        colorTicks.doLayout();
    },

    statics: {
        widgetName: 'phase_plot',
        supportsDataObjects: true,
        supportsParameterFiles: false,
        displayName: 'Phase Plot',
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
                examine = win;
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
            win.query("#weightField")[0].store = 
                ['None'].concat(obj.field_list);
            win.query("#field")[0].store = obj.field_list;
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
});
