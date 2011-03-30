function cell_finished(result, new_cell) {
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
            leaf:false, 
            expanded:true, 
            iconCls: 'pf_icon'}));
        this_pf = treePanel.root.lastChild
        Ext.each(pf.objects, function(object, obj_index) {
            this_pf.appendChild(new Ext.tree.TreeNode({text: object.name,
            leaf: true, iconCls: 'data_object'}));
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
