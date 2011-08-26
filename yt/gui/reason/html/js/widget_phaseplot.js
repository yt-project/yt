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



var WidgetPhasePlot = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "pp_" + this.id,
            title: widget_data['title'],
            iconCls: 'graph',
            autoScroll: true,
            layout:'absolute',
            closable: true,
            items: [ 
                {
                    xtype: 'panel',
                    id: 'y_ticks_' + python_varname,
                    layout: 'absolute',
                    y: 10,
                    x: 100,
                    width: 40,
                    height: 400,
                    items : [],
                    border: false,
                }, {
                    xtype: 'panel',
                    id: 'x_ticks_' + python_varname,
                    layout: 'absolute',
                    y: 410,
                    x: 140,
                    width: 400,
                    height: 40,
                    items : [],
                    border: false,
                }, {
                    xtype:'panel',
                    id: 'image_panel_' + this.id,
                    autoEl: {
                        tag: 'img',
                        id: "img_" + this.id,
                        width: 400,
                        height: 400,
                        style: 'border: 1px solid #000000',
                    },
                    x: 138,
                    y: 8,
                    width: 400,
                    height: 400,
                }, {
                    xtype:'panel',
                    id: 'colorbar_' + python_varname,
                    autoEl: {
                        tag: 'img',
                        id: "cb_" + python_varname,
                        src: "data:image/png;base64," +
                             widget_data['colorbar'],
                        width: 28,
                        height: 398,
                        style: 'border: 1px solid #000000;',
                    },
                    x: 560,
                    y: 10,
                    width: 30,
                    height: 400,
                }, {
                    xtype: 'panel',
                    id: 'ticks_' + python_varname,
                    layout: 'absolute',
                    y: 10,
                    x: 590,
                    width: 40,
                    height: 400,
                    items : [],
                    border: false,
                },{
                    xtype: 'button',
                    text: 'Upload Image',
                    x: 10,
                    y: 285,
                    width: 80,
                    tooltip: "Upload the current image to " +
                             "<a href='http://imgur.com'>imgur.com</a>",
                    handler: function(b,e) {
                        img_data = image_dom.src;
                        yt_rpc.ExtDirectREPL.upload_image(
                            {image_data:img_data,
                             caption:metadata_string},
                        function(rv) {
                            var alert_text;
                            if(rv['uploaded'] == false) {
                                alert_text = "Failure uploading image!";
                            } else {
                                alert_text = "Uploaded to " +
                                        rv['upload']['links']['imgur_page'];
                            }
                            Ext.Msg.alert('imgur.com', alert_text);
                            var record = new logging_store.recordType(
                                {record: alert_text });
                            logging_store.add(record, number_log_records++);
                        }); 
                    }
                },{
                    xtype: 'panel',
                    layout: 'vbox',
                    id: 'rhs_panel_' + python_varname,
                    width: 300,
                    height: 460,
                    x: 640, y: 10,
                    layoutConfig: {
                        align: 'stretch',
                        pack: 'start',
                    },
                    items: [
                        {
                          xtype: 'panel',
                          title: 'Plot MetaData',
                          id: 'metadata_' + python_varname,
                          style: {fontFamily: '"Inconsolata", monospace'},
                          html: 'Welcome to the Plot Window.',
                          height: 200,
                        }, {
                          xtype: 'panel',
                          title: 'Plot Editor',
                          id: 'plot_edit',
                          flex: 1,
                        }]
                }
            ]
        }
    );

    viewport.get("center-panel").activate("pp_" + this.id);
    viewport.get("center-panel").doLayout();
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("pp_" + python_varname);
    this.panel.doLayout();
    this.panel.show();
    this.image_panel = this.panel.get("image_panel_"+python_varname);
    this.ticks = this.panel.get("ticks_"+python_varname);
    var x_ticks = this.panel.get("x_ticks_"+python_varname);
    var y_ticks = this.panel.get("y_ticks_"+python_varname);
    var ticks = this.ticks;
    this.metadata_panel = this.panel.get("rhs_panel_" + python_varname).get("metadata_" + python_varname);
    var image_dom = this.image_panel.el.dom;
    var control_panel = this.panel;
    var metadata_string;
    var colorbar = this.panel.get("colorbar_"+python_varname);

    this.accept_results = function(payload) {
        this.image_panel.el.dom.src = "data:image/png;base64," + payload['image_data'];
        examine = this.metadata_panel;
        this.metadata_panel.update(payload['metadata_string']);
        metadata_string = payload['metadata_string'];
        ticks.removeAll();
        colorbar.el.dom.src = "data:image/png;base64," +
                              payload['cbar']['cmap_image'];
        Ext.each(payload['yax']['ticks'], function(tick, index) {
            examine = tick;
            y_ticks.add({xtype:'panel',
                       width: 20, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: right;' +
                              'padding-right: 5px;',
                       html: ' ' + tick[2] + ' ',
                       x:0, y: tick[0]-6});
            y_ticks.add({xtype:'panel',
                       width: 20, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:20, y: tick[0]});
        });
        y_ticks.doLayout();
        Ext.each(payload['xax']['ticks'], function(tick, index) {
            examine = tick;
            x_ticks.add({xtype:'panel',
                       width: 1, height:20,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:(400 - tick[0]) + 10, y: 0});
            x_ticks.add({xtype:'panel',
                       width: 20, height:20,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: center;',
                       html: ' ' + tick[2] + ' ',
                       x: (400 - tick[0]), y: 20});
        });
        x_ticks.doLayout();
        Ext.each(payload['cbar']['ticks'], function(tick, index) {
            ticks.add({xtype:'panel',
                       width: 10, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:0, y: tick[0]});
            ticks.add({xtype:'panel',
                       width: 30, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px;',
                       html: ' ' + tick[2] + ' ',
                       x:18, y: tick[0]-6});
        });
        ticks.doLayout();
        x_ticks.doLayout();
    }
    yt_rpc.ExtDirectREPL.execute(
        {code:python_varname + '._setup_plot()', hide:true},
        cell_finished);
}

widget_types['phase_plot'] = WidgetPhasePlot;
