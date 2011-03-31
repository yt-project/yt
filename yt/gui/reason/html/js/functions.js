function cell_finished(result, new_cell) {
    var new_log = false;
    Ext.each(result['payloads'], 
    function(payload, index) {
        if (payload['type'] == 'png_string') {
            new_cell.add(new Ext.Panel({
                autoEl:{
                    tag:'img', 
                    width:'25%',
                    src:"data:image/png;base64," + payload['image_data'],
                    id:"payload_image_" + number_images,
                    onClick: "display_image('payload_image_" + number_images + "');"
		        }
            }));
	        new_cell.doLayout();
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
            objdata: {fn: pf.filename, varname: pf.varname},
            leaf:false, 
            expanded:true, 
            iconCls: 'pf_icon'}));
        this_pf = treePanel.root.lastChild
        Ext.each(pf.objects, function(obj, obj_index) {
            this_pf.appendChild(new Ext.tree.TreeNode(
                {text: obj.name,
                 leaf: true,
                 iconCls: 'data_obj',
                 objdata: {varname: obj.varname},
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
        width:240,
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
                id: 'x_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Y',
                id: 'y_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
                xtype:'textfield',
                fieldLabel: 'Center Z',
                id: 'z_center',
                value: '0.5',
                width: 90,
                allowBlank:false,
            },{
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
                store:['Density','Temperature','X Velocity','Y Velocity','Z Velocity'],
                width: 90,
                allowBlank:false,
            }],
            buttons: [
                {
                    text: 'Slice',
                    handler: function(b, e){Ext.Msg.alert('Slicing','Slicing it up!')}
                },{
                    text: 'Cancel',
                    handler: function(b, e){Ext.Msg.alert('Cancelled','Slice cancelled.')}
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
        {code: fcall}, handle_payload);
}


function getProjectionHandler(node){
function projectionHandler(item,pressed){
    var win = new Ext.Window({
        layout:'fit',
        width:240,
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
                store:['Density','Temperature','X Velocity','Y Velocity','Z Velocity'],
                width: 120,
                allowBlank:false,
            },{
                xtype:'combo',
                fieldLabel: 'Weight Field',
                id: 'weightField',
                store:['None','Density','Temperature','X Velocity','Y Velocity','Z Velocity'],
                width: 120,
                allowBlank:false,
            }],
            buttons: [
                {
                    text: 'Project',
                    handler: function(b, e){Ext.Msg.alert('Projection','Projecting!')}
                },{
                    text: 'Cancel',
                    handler: function(b, e){Ext.Msg.alert('Cancelled','Projection cancelled.')}
                }
            ]
        }]
    });
    win.show(this);
}
return projectionHandler;
}
=======
>>>>>>> other
