/**********************************************************************
The Plot Window Widget

Author: Cameron Hummels <chummels@gmail.com>
Affiliation: Columbia
Author: Jeffrey S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
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



var WidgetPlotWindow = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;
    this.print_python = function(b, e) {
        yt_rpc.ExtDirectREPL.execute(
            {code:'print "' + python_varname + '"',
             hide:true},
            function(f, a) {alert(a.result['output']);}
        );
    }

    this.widget_keys = new Ext.KeyMap(document, [
        {key: 'z', fn: function(){control_panel.get("zoom10x").handler();}}
    ]);
    var widget_keys = this.widget_keys;
    widget_keys.disable();
    widget_keys.varname = python_varname;

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "pw_" + this.id,
            title: "Plot Window",
            iconCls: 'graph',
            autoScroll: true,
            layout:'absolute',
            closable: true,
            listeners: {activate: function(p){
                                widget_keys.enable();
                            },
                        deactivate: function(p){
                                widget_keys.disable();
                            }
                        },
            items: [ 
                {
                    xtype:'panel',
                    id: 'image_panel_' + this.id,
                    autoEl: {
                        tag: 'img',
                        id: "img_" + this.id,
                        width: 400,
                        height: 400,
                    },
                    x: 100,
                    y: 10,
                    width: 400,
                    height: 400,
                }, 
                /* the single buttons for 10% pan*/
                {
                    xtype:'button',
                    iconCls: 'singleuparrow',
                    //text: 'North',
                    x: 40,
                    y: 10,
                    handler: function(b,e) {
                        cc = python_varname + '.pan_rel((0.0, -0.1))'
                        yt_rpc.ExtDirectREPL.execute(
                        {code:cc, hide:true}, cell_finished); 
                    }
                }, {
                    xtype:'button',
                    iconCls: 'singlerightarrow',
                    //text:'East',
                    x : 60,
                    y : 30,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.1, 0.0))',
                             hide:true},
                        cell_finished); 
                    }
                }, {
                    xtype:'button',
                    iconCls: 'singledownarrow',
                    //text: 'South',
                    x: 40,
                    y: 50,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.0, 0.1))',
                             hide:true},
                        cell_finished); 
                    }
                }, {
                    xtype: 'button',
                    iconCls: 'singleleftarrow',
                    //text: 'West',
                    x: 20,
                    y: 30,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((-0.1, 0.0))',
                             hide:true},
                        cell_finished); 
                    }
                }, 
                /* the double buttons for 50% pan*/
                {
                    xtype:'button',
                    iconCls: 'doubleuparrow',
                    //text: 'North',
                    x: 40,
                    y: 80,
                    handler: function(b,e) {
                        cc = python_varname + '.pan_rel((0.0, -0.5))'
                        yt_rpc.ExtDirectREPL.execute(
                        {code:cc, hide:true}, cell_finished); 
                    }
                }, {
                    xtype:'button',
                    iconCls: 'doublerightarrow',
                    //text:'East',
                    x : 60,
                    y : 100,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.5, 0.0))',
                             hide:true},
                        cell_finished); 
                    }
                }, {
                    xtype:'button',
                    iconCls: 'doubledownarrow',
                    //text: 'South',
                    x: 40,
                    y: 120,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.0, 0.5))',
                             hide:true},
                        cell_finished); 
                    }
                }, {
                    xtype: 'button',
                    iconCls: 'doubleleftarrow',
                    //text: 'West',
                    x: 20,
                    y: 100,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((-0.5, 0.0))',
                             hide:true},
                        cell_finished); 
                    }
                },
                /* Now the zoom buttons */
                {
                    xtype: 'button',
                    text: 'Zoom In 10x',
                    id: "zoom10x",
                    x: 10,
                    y: 160,
                    width: 80,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.zoom(10.0)',
                             hide:true},
                        cell_finished); 
                    }
                },{
                    xtype: 'button',
                    text: 'Zoom In 2x',
                    x: 10,
                    y: 185,
                    width: 80,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.zoom(2.0)',
                             hide:true},
                        cell_finished); 
                    }
                },{
                    xtype: 'button',
                    text: 'Zoom Out 2x',
                    x: 10,
                    y: 210,
                    width: 80,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.zoom(0.5)',
                             hide:true},
                        cell_finished); 
                    }
                },{
                    xtype: 'button',
                    text: 'Zoom Out 10x',
                    x: 10,
                    y: 235,
                    width: 80,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.zoom(0.1)',
                             hide:true},
                        cell_finished); 
                    }
                },{
                    xtype: 'button',
                    text: 'Upload Image',
                    x: 10,
                    y: 285,
                    width: 80,
                    handler: function(b,e) {
                        img_data = image_dom.src;
                        yt_rpc.ExtDirectREPL.upload_image(
                            {image_data:img_data},
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
                    xtype: 'combo',
                    text: 'Field',
                    x: 10,
                    y: 315,
                    width: 80,
                    store:widget_data['fields'],
                    value:widget_data['initial_field'],
                    editable: false,
                    triggerAction: 'all',
                    listeners: {change: function(form, newValue, oldValue) {
                        alert(newValue);
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.set_current_field("' +
                                newValue + '")', hide:true},
                            cell_finished);
                    }}
                }
            ]
        }
    );

    viewport.get("center-panel").activate("pw_" + this.id);
    viewport.get("center-panel").doLayout();
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("pw_" + python_varname);
    this.panel.doLayout();
    this.panel.show();
    this.image_panel = this.panel.get("image_panel_"+python_varname);
    var image_dom = this.image_panel.el.dom;
    var control_panel = this.panel;

    this.accept_results = function(payload) {
        this.image_panel.el.dom.src = "data:image/png;base64," + payload['image_data'];
    }

    yt_rpc.ExtDirectREPL.execute(
        {code:python_varname + '.zoom(1.0)', hide:true},
        cell_finished);
}

widget_types['plot_window'] = WidgetPlotWindow;
