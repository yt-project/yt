/**********************************************************************
The main GUI facility for Reason

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


var viewport;
var widget_types = {}
var widget_list = {}

var examine;
var number_log_records = 0;
var number_images = 0;

var res;
var cell_count = 0;

var handle_result = function(f, a) {
    var input_line = repl_input.get("input_line")
    if (a.result == null) {
        formatted_input = input_line.getValue();
        text = "ERROR";
        var cell = new_cell(formatted_input, text);
        OutputContainer.add(cell);
        notebook.doLayout();
        input_line.setValue("");
        return;
    }
    cell_finished(a.result);
}

    var handle_payload = function(pp) {
        cell_finished(pp);
    }

    var repl_input = new Ext.FormPanel({
        title: 'YT Input',
        url: 'push',
        flex: 0.2,
        layout: 'fit',
        padding: 5,
        items: [{
            id: 'input_line',
            xtype: 'textarea',
            width: '100%',
            autoScroll: true,
            name: 'line',
            allowBlank: 'True',
            bodyStyle: 'font-family: "monospace";',
            listeners: {
                specialkey: function(f, e){
                    if (e.getKey() == e.ENTER) {
                        cell_sent();
                        yt_rpc.ExtDirectREPL.execute({
                            code:repl_input.get('input_line').getValue()},
                        handle_result);
	                }
	            },
                afterrender: function(f, e){
                    //var input_line_drop_target_el = repl_input.get("input_line").el.dom;
                    var input_line_drop_target_el = repl_input.body.dom;

                    var input_line_drop_target = new Ext.dd.DropTarget(input_line_drop_target_el, {
                        ddGroup     : 'pfDDgroup',
                        notifyEnter : function(ddSource, e, data) {
                            repl_input.body.stopFx();
                            repl_input.body.highlight();
                        },
                        notifyDrop  : function(ddSource, e, data){

                            var varname = data.node.attributes.objdata.varname;
                            /* There is possibly a better way to do this, where it's also
                               inserted correctly. */
                            var line = repl_input.get("input_line");
                            line.setValue(line.getValue() + varname);
                            line.focus();
                            return(true);
                        }
                    });
                },
            },
        },],
    });

    var NorthButton = new Ext.Button({
        text : 'North',
	    pageX : 205,
        pageY : 10
//        handler: function(b, e) { window.open("session.py", "_top"); }
    });

    var EastButton = new Ext.Button({
        text:'East',
        pageX : 410,
        pageY : 205
    });

    var SouthButton = new Ext.Button({
        text:'South',
        pageX : 205,
        pageY : 410
    });

    var WestButton = new Ext.Button({
        text:'West',
        pageX : 10,
        pageY : 205
    });

    var OutputContainer = new Ext.Panel({
        title: 'YT Output',
        id: 'output_container',
        autoScroll: true,
        flex: 0.8,
        items: []
    });

    var PlotPanel = new Ext.Panel({
        title: 'Plot Window 1',
        iconCls: 'graph',
        autoScroll: true,
        layout:'absolute',
        items: [ 
            NorthButton,
            EastButton,
            SouthButton,
            WestButton
        ]
    });

    var examine;
    var notebook;


    var treePanel = new Ext.tree.TreePanel({
        iconCls: 'nav',
        id: 'tree-panel',
        layout: 'anchor',
        region:'west',
        split: true,
        anchor: '100% -35',
        minSize: 150,
        autoScroll: true,
        rootVisible: false,
        ddGroup: 'pfDDgroup',
        enableDD: true,
        root:new Ext.tree.TreeNode({
            expanded:true,
            leaf:false,
            text:''
        }),
        listeners: {
                render: {
                    fn: function() {
                        Ext.getBody().on("contextmenu", Ext.emptyFn,
                            null, {preventDefault: true});
                    }
                },
                contextmenu: {
                    fn: function(node, event){
                        if (node.attributes.objdata.type == 'obj') return;
                        var rightClickMenu = new Ext.menu.Menu({
                            items: [
                                {
                                    text: 'Open slice',
                                    handler: getSliceHandler(node),
                                }, {
                                    text: 'Open projection',
                                    handler: getProjectionHandler(node),
                                }
                            ]
                        });
                        rightClickMenu.showAt(event.xy);
                    }
                }
          }
    });

    var ButtonGroupPanel = new Ext.Panel({
        layout: 'anchor',
        ButtonAlign: 'center',
        collapsible: false,
        renderTo: document.body,
        tbar: [{
            xtype: 'buttongroup',
            columns: 7,
            items: [{
                text: 'Download',
                layout:'anchor',
                anchor: '100% 25%',
                handler: function(b, e) { 
                    window.open("session.py", "_top"); 
                    var record = new logging_store.recordType({
                        record: 'Saved session locally.'});
                    logging_store.add(record, number_log_records++);
	            }
            },{
                xtype: 'tbseparator'
            },{
                text: 'Save',
                layout:'anchor',
	            anchor: '100% 50%',
	            handler: function (b,e) { 
                    Ext.Msg.prompt("We have important work to do.", 
                    "Enter filename.", 
                    function(btn, text) {
                        if (btn == 'ok'){
                            yt_rpc.ExtDirectREPL.save_session({filename:text}, 
                            function(f, a) {
                                if (a.result['status'] == 'SUCCESS') {
                                    var alert_text = 'Saved session to ' + 
                                    a.result['filename']
                                    Ext.Msg.alert('Success!', alert_text);
                                    var record = new logging_store.recordType(
                                        {record: alert_text });
                                    logging_store.add(record, number_log_records++);
							    } else {
							        Ext.Msg.alert('Always naysaying!',
                                        'Failed to save to ' + 
                                        a.result['filename'] + 
                                        '<br>Error: ' + 
                                        a.result['error']);
                                }
                            });
                        }
                    });
                }
            },{
                xtype: 'tbseparator'
            },{
                text: 'Paste',
                layout:'anchor',
                anchor: '100% 75%',
                handler: function (b,e) { 
                    yt_rpc.ExtDirectREPL.paste_session({}, function(f, a) {
                        if (a.result['status'] == 'SUCCESS') {
                            var alert_text = 'Pasted session to:<br>' + 
                            a.result['site']
                            var alert_text_rec = 'Pasted session to: ' + 
                            a.result['site']
                            Ext.Msg.alert('Pastebin', alert_text);
                            var record = new logging_store.recordType(
                                {record: alert_text_rec });
                            logging_store.add(record, number_log_records++);
                        }
                    }); 
                }
            },{
                xtype: 'tbseparator'
            },{
                text: 'Help',
                layout:'anchor',
                anchor: '100% 100%',
                handler: function (b,e) { 
                        window.open("help.html", "_new");
                }
            }]
        }]
    });

    var status_panel;
    var logging_store = new Ext.data.Store({
        fields: [{name:'record'}],
        reader: new Ext.data.ArrayReader({}, [{name: 'record'}]),
    });

    var heartbeep_request = false;
    var task_runner = new Ext.util.TaskRunner();


    Ext.onReady(function(){
       Ext.BLANK_IMAGE_URL = 'resources/resources/images/default/s.gif';

    // NOTE: This is an example showing simple state management. During development,
    // it is generally best to disable state management as dynamically-generated ids
    // can change across page loads, leading to unpredictable results.  The developer
    // should ensure that stable state ids are set for stateful components in real apps.
    // it's a cold day for pontooning.
        Ext.state.Manager.setProvider(new Ext.state.CookieProvider());

    // Go ahead and create the TreePanel now so that we can use it below
        viewport = new Ext.Viewport({
            layout: 'border',
            items: [
		    // lazily created panel (xtype:'panel' is default)
                {
                    xtype: 'grid',
                    store: logging_store,
                    defaults: { width: 800 },
                    columns: [ {id:'record', 
                        sortable: false,
                        width:800} ],
                    autofill: true,
                    region: 'south',
                    id: "status-region",
                    cls: "status-logger",
                    split: true,
                    height: 100,
                    maxSize: 200,
                    collapsible: true,
                    title: 'Status',
                    margins: '0 0 0 0',
                }, {
                    region: 'west',
                    id: 'west-panel', // see Ext.getCmp() below
                    title: 'Data Objects',
                    split: true,
                    width: 200,
                    minSize: 175,
                    maxSize: 400,
                    collapsible: true,
                    margins: '0 0 0 5',
                    layout: {
                        type: 'anchor',
                    },
                    items: [
                        treePanel,
                        ButtonGroupPanel
                    ]
		  // in this instance the TabPanel is not wrapped by another panel
		  // since no title is needed, this Panel is added directly
		  // as a Container
                },{
                    xtype: 'tabpanel',
                    region: 'center', 
                    id: 'center-panel',
                    deferredRender: false,
                    activeTab: 0,     
                    items: [
                        {
                            title: 'YT',
                            id: 'notebook',
                            layout: 'vbox',
                            layoutConfig: {align:'stretch'},
                            closable: false,
                            autoScroll: false,
                            iconCls: 'console',
                            items: [repl_input, OutputContainer]
                        }, 
//                        PlotPanel
                    ]
                }
            ]
        });

// get a reference to the HTML element with id "hideit" and add a click listener to it 
    console.log('Mitchell!\nPardon me! Mitchell!')
    Ext.get("hideit").on('click', function(){
// get a reference to the Panel that was created with id = 'west-panel' 
	    var w = Ext.getCmp('west-panel');
// expand or collapse that Panel based on its collapsed property state
// need to make room for six sour cream burritos
        w.collapsed ? w.expand() : w.collapse();
    });
    
    notebook = viewport.get("center-panel").get("notebook");
    status_panel = viewport.get("status-region").get("status-div");
    
    var record = new logging_store.recordType(
        {record: 'Welcome to yt.'});
    logging_store.add(record, number_log_records++);

    var record = new logging_store.recordType(
        {record: 'After entering a line of code in the YT Input field, press shift-enter to evaluate.' });
    logging_store.add(record, number_log_records++);

    var record = new logging_store.recordType(
        {record: '4d3d3d3 engaged.' });
    logging_store.add(record, number_log_records++);

    if (!Ext.state.Manager.get("reason_welcomed", false)) {
        Ext.MessageBox.alert("Reason v0.5",
        "Welcome to Reason.  <br>Treat the 'YT Input' field as a YT/python intepreter.<br>Press shift-enter to evaluate.",
        function(b,e){ repl_input.get("input_line").focus(); });
        Ext.state.Manager.set("reason_welcomed", true);
    } else { 
        repl_input.get("input_line").focus(); }

    /* Set up the heartbeep */
    var heartbeep = {
    run:
      function(){ if (heartbeep_request == true) return; 
        yt_rpc.ExtDirectREPL.heartbeep(
            {}, function(f, a) {
            heartbeep_request = false;
            if (f['alive'] == true) { 
            }})},
    interval: 5000};

    //task_runner.start(heartbeep);
                         
    });
