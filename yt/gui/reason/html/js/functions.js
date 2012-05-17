/**********************************************************************
Functions for Reason

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

function cell_finished(result) {
    var new_log = false;
    var cell_resulted = false;
    var cell;
    Ext.each(result, 
    function(payload, index) {
        if (payload['type'] == 'shutdown') {
            reason.task_runner.stop(heartbeat);
            reason.heartbeat_request = true;
            return;
        } else if (payload['type'] == 'cell_results') {
            text = "<pre>"+payload['output']+"</pre>";
            formatted_input = payload['input']
            cell = new_cell(formatted_input, text, payload['raw_input']);
            OutputContainer.add(cell);
            OutputContainer.doLayout();
            notebook.doLayout();
            if (repl_input.locked == true) {
                /* Assume only one locking level */
                repl_input.locked = false;
            } else {
                repl_input.get("input_line").setValue("");
            }
            if (OutputContainer.items.length > 1) {
                OutputContainer.body.dom.scrollTop = 
                OutputContainer.body.dom.scrollHeight -
                cell.body.dom.scrollHeight - 20;
            }
            cell_resulted = true;
        } else if (payload['type'] == 'png_string') {
            OutputContainer.add(new Ext.Panel({
                autoEl:{
                    tag:'img', 
                    width:'25%',
                    src:"data:image/png;base64," + payload['image_data'],
                    id:"payload_image_" + number_images,
                    onClick: "display_image('payload_image_" + number_images + "');"
		        }
            }));
	        OutputContainer.doLayout();
	        number_images++;
        } else if (payload['type'] == 'cell_contents') {
	        var input_line = repl_input.get("input_line");
	        input_line.setValue(payload['value']);
            repl_input.locked = true;
        } else if (payload['type'] == 'log_entry') {
            reason.log(payload['log_entry']);
            new_log = true;
        } else if (payload['type'] == 'widget') {
            var widget_type = payload['widget_type'];
            var widget = new widget_types[widget_type](payload['varname'],
                                                       payload['data']);
            widget_list[widget.id] = widget;
            /*
               Sometimes instantiating a widget adds some objects ...
               Plus, often when creating a widget we disable the 
               entry of data and whatnot. 
            */
            cell_resulted = true;
        } else if (payload['type'] == 'widget_payload') {
            var widget = widget_list[payload['widget_id']];
            widget.accept_results(payload);
        } else {
            alert("Didn't know how to process " + payload['type']);
        }
    });
    if (new_log == true){
        reason.log_scroll()
    }
    if (cell_resulted == true) {
        enable_input();
    }
}

function display_image(image_id) {
    var image = Ext.get(image_id);
    var src = image.dom.src;
    var virtualdom = '<html><title>Image Viewer</title><body><img src="' 
        + src + '"/></body></html>',
    prev = window.open('', 'image_viewer');
    prev.document.open();
    prev.document.write(virtualdom);
    prev.document.close();
}

// Create a tree in the left panel with the pfs and their objects.
function fill_tree(my_pfs) {
    treePanel.root.removeAll();
    Ext.each(my_pfs, function(pf, index) {
        treePanel.root.appendChild(new Ext.tree.TreeNode({
            text: pf.name,
            objdata: {fn: pf.filename, varname: pf.varname, type: 'pf',
                      field_list: pf.field_list},
            leaf:false, 
            expanded:true, 
            iconCls: 'pf_icon'}));
        this_pf = treePanel.root.lastChild
        Ext.each(pf.objects, function(obj, obj_index) {
            this_pf.appendChild(new Ext.tree.TreeNode(
                {text: obj.name,
                 leaf: true,
                 iconCls: 'data_obj',
                 objdata: {varname: obj.varname, type: 'obj',
                           pfdata: this_pf.attributes.objdata},
                 }));
        });
    });
}

function new_cell(input, result, raw_input) {
    var name = "cell_" + cell_count;
    var CellPanel = new Ext.Panel(
        { 
            id: name, 
            //title: "Cell " + cell_count,
            items: [
                { xtype:'panel',
                  layout: 'hbox',
                  id:name+"_input",
                  items: [
                    { xtype:'panel',
                      html:input,
                      flex:1,
                      boxMinHeight: 40,
                    },
                    { xtype: 'button',
                      width: 24,
                      height: 24,
                      iconCls: 'upload',
                      tooltip: 'Upload to Pastebin',
                      listeners: {
                          click: function(f, e) {
                            yt_rpc.ExtDirectREPL.paste_text({to_paste:raw_input},
                              function(f, a) {
                                if (a.result['status'] == 'SUCCESS') {
                                    var alert_text = 'Pasted cell to:<br>' + 
                                    a.result['site']
                                    var alert_text_rec = 'Pasted cell to: ' + 
                                    a.result['site']
                                    Ext.Msg.alert('Pastebin', alert_text);
                                    var record = new logging_store.recordType(
                                        {record: alert_text_rec });
                                    logging_store.add(record, number_log_records++);
                              }
                            });
                          }
                        }
                    },
                    { xtype: 'button',
                      width: 24,
                      height: 24,
                      iconCls: 'doubleuparrow',
                      tooltip: 'Copy into current cell',
                      listeners: {
                          click: function(f, e) {
                            repl_input.get('input_line').setValue(raw_input);
                          }
                      },
                    },
                  ],
                },
                { xtype:'panel',
                  layout: 'hbox',
                  items: [
                    { xtype:'panel',
                      id:name+"_result",
                      autoScroll:true,
                      flex: 1,
                      html:result,
                      boxMinHeight: 40,
                    },
                  ],
                },
            ]
        }
    );
    cell_count++;
    return CellPanel;
}

function getGridViewerHandler(node){
function gridViewerHandler(item, pressed){
    yt_rpc.ExtDirectREPL.create_grid_viewer(
        {pfname:node.attributes.objdata.varname},
        handle_result);
}
return gridViewerHandler;
}

function getGridDataViewerHandler(node){
function gridDataViewerHandler(item, pressed){
    yt_rpc.ExtDirectREPL.create_grid_dataview(
        {pfname:node.attributes.objdata.varname},
        handle_result);
}
return gridDataViewerHandler;
}

function getStreamlineViewerHandler(node){
function streamlineViewerHandler(item, pressed){
    yt_rpc.ExtDirectREPL.create_streamline_viewer(
        {pfname:node.attributes.objdata.varname},
        handle_result);
}
return streamlineViewerHandler;
}

function getIsocontourViewerHandler(node){
function isocontourViewerHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:320,
        height:250,
        modal:true,
        resizable:false,
        draggable:false,
        border:false,
        title:'Isocontour Extraction in' + node,
        items: [{
            xtype: 'form', // FormPanel
            labelWidth:80,
            frame:true,
            items: [{
                xtype:'combo',
                fieldLabel: 'Field',
                id: 'field',
                store:node.attributes.objdata.field_list,
                width: 200,
                allowBlank:false,
                value: 'Density',
                triggerAction: 'all',
            },{
                xtype:'combo',
                fieldLabel: 'Sampling Field',
                id: 'extract_field',
                store:node.attributes.objdata.field_list,
                width: 200,
                allowBlank:false,
                value: 'Temperature',
                triggerAction: 'all',
            },{
                xtype:'textfield',
                fieldLabel: 'Value',
                id: 'value',
                value: '1e-25',
                width: 90,
                allowBlank:false,
            }],
            buttons: [
                {
                    text: 'Extract',
                    handler: function(b, e){
                        var field = Ext.get("field").getValue();
                        var value = Ext.get("value").getValue();
                        var sampling_field = Ext.get("extract_field").getValue();
                        yt_rpc.ExtDirectREPL.create_isocontours({
                            pfname:node.attributes.objdata.varname,
                            field:field, value:value,
                            sampling_field:sampling_field},
                          handle_result);
                        disable_input();
                        win.close();
                    }
                },{
                    text: 'Cancel',
                    handler: function(b, e){
                        win.close();
                    }
                }
            ]
        }]
    });
    win.show(this);
}
return isocontourViewerHandler;
}

function getSliceHandler(node){
function sliceHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:320,
        height:250,
        modal:true,
        resizable:false,
        draggable:false,
        border:false,
        title:'Slice Details for ' + node,
        items: [{
            xtype: 'form', // FormPanel
            labelWidth:80,
            frame:true,
            items: [{
                xtype:'textfield',
                fieldLabel: 'Center X',
                id: 'slice_x_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Y',
                id: 'slice_y_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Z',
                id: 'slice_z_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'combo',
                fieldLabel: 'Axis',
                id: 'slice_axis',
                store:['X','Y','Z'],
                width: 90,
                allowBlank:false,
                value: 'X',
                triggerAction: 'all',
            },{
                xtype:'checkbox',
                fieldLabel: 'Center on Max',
                id: 'max_dens',
                width: 90,
                allowBlank:false,
                handler: function(checkbox, checked) {
                    if (checked == true) {
                        this.ownerCt.get("slice_x_center").disable();
                        this.ownerCt.get("slice_y_center").disable();
                        this.ownerCt.get("slice_z_center").disable();
                    } else {
                        this.ownerCt.get("slice_x_center").enable();
                        this.ownerCt.get("slice_y_center").enable();
                        this.ownerCt.get("slice_z_center").enable();
                    }
                }
            },{
                xtype:'combo',
                fieldLabel: 'Field',
                id: 'slice_field',
                store:node.attributes.objdata.field_list,
                width: 200,
                allowBlank:false,
                value: 'Density',
                triggerAction: 'all',
            }],
            buttons: [
                {
                    text: 'Slice',
                    handler: function(b, e){
                        var center = [Ext.get("slice_x_center").getValue(),
                                      Ext.get("slice_y_center").getValue(),
                                      Ext.get("slice_z_center").getValue()];
                        var axis = Ext.get("slice_axis").getValue();
                        var field = Ext.get("slice_field").getValue();
                        var onmax = Ext.get("max_dens").getValue();
                        yt_rpc.ExtDirectREPL.create_slice({
                            pfname:node.attributes.objdata.varname,
                            center: center, axis:axis, field:field, onmax:onmax},
                          handle_result);
                        disable_input();
                        win.close();
                    }
                },{
                    text: 'Cancel',
                    handler: function(b, e){
                        win.close();

                    }
                }
            ]
        }]
    });
    win.show(this);
}
return sliceHandler;
}

function widget_call(varname, method) {
    var fcall = varname + "." + method;
    yt_rpc.ExtDirectREPL.execute(
        {code: fcall}, cell_finished);
}

function getPhasePlotHandler(node){
function phasePlotHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:370,
        height:220,
        modal:true,
        resizable:false,
        draggable:false,
        border:false,
        title:'Phase Plot Details for ' + node,
        items: [{
            xtype: 'form', // FormPanel
            labelWidth:80,
            frame:true,
            items: [ {
                xtype:'combo',
                fieldLabel: 'X Field',
                id: 'x_field',
                store:node.attributes.objdata.pfdata.field_list,
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'Density'
            },{
                xtype:'combo',
                fieldLabel: 'Y Field',
                id: 'y_field',
                store:node.attributes.objdata.pfdata.field_list,
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'Temperature'
            },{
                xtype:'combo',
                fieldLabel: 'Z Field',
                id: 'z_field',
                store:node.attributes.objdata.pfdata.field_list,
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'CellMassMsun'
            },{
                xtype:'combo',
                fieldLabel: 'Weight Field',
                id: 'weight',
                store:['None'].concat(node.attributes.objdata.pfdata.field_list),
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'None'
            }],
            buttons: [
                {
                    text: 'Calculate',
                    handler: function(b, e){
                        var x_field = Ext.get("x_field").getValue();
                        var y_field = Ext.get("y_field").getValue();
                        var z_field = Ext.get("z_field").getValue();
                        var weight = Ext.get("weight").getValue();
                        yt_rpc.ExtDirectREPL.create_phase({
                                objname: node.attributes.objdata.varname,
                                /* Mirror image varnames ... */
                                field_x: x_field,
                                field_y: y_field,
                                field_z: z_field,
                                weight: weight,
                                },
                              handle_result);
                        disable_input();
                        win.close();
                    }
                },{
                    text: 'Cancel',
                    handler: function(b, e){win.close()}
                }
            ]
        }]
    });
    win.show(this);
}
return phasePlotHandler;
}

function getProjectionHandler(node){
function projectionHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:370,
        height:220,
        modal:true,
        resizable:false,
        draggable:false,
        border:false,
        title:'Projection Details for ' + node,
        items: [{
            xtype: 'form', // FormPanel
            labelWidth:80,
            frame:true,
            items: [{
                xtype:'combo',
                fieldLabel: 'Axis',
                id: 'axis',
                store:['X','Y','Z'],
                width: 90,
                allowBlank:false,
                triggerAction: 'all',
                value: 'X',
            },{
                xtype:'checkbox',
                fieldLabel: 'Center on Max',
                id: 'max_dens',
                width: 90,
                allowBlank:false,
                /* No handler, because no center */
            },{
                xtype:'combo',
                fieldLabel: 'Field',
                id: 'field',
                store:node.attributes.objdata.field_list,
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'Density'
            },{
                xtype:'combo',
                fieldLabel: 'Weight Field',
                id: 'weightField',
                store:['None'].concat(node.attributes.objdata.field_list),
                width: 230,
                allowBlank:false,
                triggerAction: 'all',
                value: 'None'
            }],
            buttons: [
                {
                    text: 'Project',
                    handler: function(b, e){
                        var axis = Ext.get("axis").getValue();
                        var field = Ext.get("field").getValue();
                        var weight = Ext.get("weightField").getValue();
                        var onmax = Ext.get("max_dens").getValue();
                        yt_rpc.ExtDirectREPL.create_proj({
                                pfname: node.attributes.objdata.varname,
                                axis: axis, field: field, weight: weight,
                                onmax: onmax},
                              handle_result);
                        disable_input();
                        win.close();
                    }
                },{
                    text: 'Cancel',
                    handler: function(b, e){win.close()}
                }
            ]
        }]
    });
    win.show(this);
}
return projectionHandler;
}

function getSphereCreator(node){
function sphereCreator(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:320,
        height:250,
        modal:true,
        resizable:false,
        draggable:false,
        border:false,
        title:'Sphere Creator ' + node,
        items: [{
            xtype: 'form', // FormPanel
            labelWidth:80,
            frame:true,
            items: [{
                xtype:'textfield',
                fieldLabel: 'Center X',
                id: 'slice_x_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Y',
                id: 'slice_y_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Z',
                id: 'slice_z_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Radius',
                id: 'radius_value',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'combo',
                fieldLabel: 'Unit',
                id: 'radius_unit',
                store:['unitary', '1', 'mpc', 'kpc', 'pc', 'au', 'rsun', 'cm'],
                width: 90,
                allowBlank:false,
                value: 'Unitary',
                triggerAction: 'all',
            },{
                xtype:'checkbox',
                fieldLabel: 'Center on Max',
                id: 'max_dens',
                width: 90,
                allowBlank:false,
                handler: function(checkbox, checked) {
                    if (checked == true) {
                        this.ownerCt.get("slice_x_center").disable();
                        this.ownerCt.get("slice_y_center").disable();
                        this.ownerCt.get("slice_z_center").disable();
                    } else {
                        this.ownerCt.get("slice_x_center").enable();
                        this.ownerCt.get("slice_y_center").enable();
                        this.ownerCt.get("slice_z_center").enable();
                    }
                }
            }],
            buttons: [
                {
                    text: 'Slice',
                    handler: function(b, e){
                        var center = [Ext.get("slice_x_center").getValue(),
                                      Ext.get("slice_y_center").getValue(),
                                      Ext.get("slice_z_center").getValue()];
                        var onmax = Ext.get("max_dens").getValue();
                        var radius = [Ext.get("radius_value").getValue(),
                                      Ext.get("radius_unit").getValue()]
                        objargs = {radius: radius}
                        if (onmax == true) {
                            objargs['center'] = 'max';
                        } else {
                            objargs['center'] = center;
                        }
                        yt_rpc.ExtDirectREPL.object_creator({
                            pfname:node.attributes.objdata.varname,
                            objtype:'sphere', objargs:objargs},
                          handle_result);
                        disable_input();
                        win.close();
                    }
                },{
                    text: 'Cancel',
                    handler: function(b, e){
                        win.close();

                    }
                }
            ]
        }]
    });
    win.show(this);
}
return sphereCreator;
}
