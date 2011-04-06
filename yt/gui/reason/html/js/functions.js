/**********************************************************************
Functions for Reason

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

function cell_finished(result) {
    var new_log = false;
    Ext.each(result, 
    function(payload, index) {
        if (payload['type'] == 'cell_results') {
            text = "<pre>"+payload['output']+"</pre>";
            formatted_input = payload['input']
            var cell = new_cell(formatted_input, text);
            OutputContainer.add(cell);
            OutputContainer.doLayout();
            notebook.doLayout();
            repl_input.get("input_line").setValue("");
            if (OutputContainer.items.length > 1) {
                examine = cell;
                OutputContainer.body.dom.scrollTop = 
                OutputContainer.body.dom.scrollHeight -
                cell.body.dom.scrollHeight - 20;
            }
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
        } else if (payload['type'] == 'log_entry') {
	        var record = new logging_store.recordType(
		        {record: payload['log_entry'] });
	        logging_store.add(record, number_log_records++);
            new_log = true;
        } else if (payload['type'] == 'widget') {
            var widget_type = payload['widget_type'];
            var widget = new widget_types[widget_type](payload['varname']);
            widget_list[widget.id] = widget;
        } else if (payload['type'] == 'widget_payload') {
            var widget = widget_list[payload['widget_id']];
            widget.accept_results(payload);
        }
    });
    yt_rpc.ExtDirectParameterFileList.get_list_of_pfs({}, fill_tree);
    if (new_log == true){
        viewport.get("status-region").getView().focusRow(number_log_records-1);
    }
    repl_input.body.removeClass("cell_waiting");
    repl_input.get('input_line').setReadOnly(false);
    repl_input.get("input_line").focus();
}

function cell_sent() {
    repl_input.get('input_line').setReadOnly(true);
    repl_input.body.addClass("cell_waiting");
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
                 objdata: {varname: obj.varname, type: 'obj'},
                 }));
        });
    });
}

function new_cell(input, result) {
    var name = "cell_" + cell_count;
    var CellPanel = new Ext.Panel(
        { 
            id: name, 
            //title: "Cell " + cell_count,
            items: [
                new Ext.Panel({
                    id:name+"_input",
                    html:input,
                }),
                new Ext.Panel({
                    id:name+"_result",
                    autoScroll:true,
                    width: "100%",
                    html:result,
                })
            ]
        }
    );
    cell_count++;
    return CellPanel;
}

function getSliceHandler(node){
function sliceHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:320,
        height:200,
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
            },{
                xtype:'combo',
                fieldLabel: 'Field',
                id: 'slice_field',
                store:node.attributes.objdata.field_list,
                width: 200,
                allowBlank:false,
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
                        yt_rpc.ExtDirectREPL.create_slice({
                            pfname:node.attributes.objdata.varname,
                            center: center, axis:axis, field:field},
                          handle_result);
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


function getProjectionHandler(node){
function projectionHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:370,
        height:170,
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
            },{
                xtype:'combo',
                fieldLabel: 'Field',
                id: 'field',
                store:node.attributes.objdata.field_list,
                width: 230,
                allowBlank:false,
            },{
                xtype:'combo',
                fieldLabel: 'Weight Field',
                id: 'weightField',
                store:['None'].concat(node.attributes.objdata.field_list),
                width: 230,
                allowBlank:false,
            }],
            buttons: [
                {
                    text: 'Project',
                    handler: function(b, e){
                        var axis = Ext.get("axis").getValue();
                        var field = Ext.get("field").getValue();
                        var weight = Ext.get("weightField").getValue();
                        yt_rpc.ExtDirectREPL.create_proj({
                                pfname: node.attributes.objdata.varname,
                                axis: axis, field: field, weight: weight},
                              handle_result);
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
