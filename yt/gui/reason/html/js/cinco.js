var viewport;

var examine;
var number_log_records = 0;
var number_images = 0;

var res;
var cell_count = 0;

var handle_result = function(f, a) {
  var input_line = repl_input.get("input_line")
  if (a.result == null) {
    text = "ERROR";
    formatted_input = input_line.getValue();
  } else {
    //text = a.result['output'].replace(/\n/g,"<br/>");
    text = "<pre>"+a.result['output']+"</pre>";
    formatted_input = a.result['input']
  }
  var cell = new_cell(formatted_input, text);
  OutputContainer.add(cell);
  notebook.doLayout();
  input_line.setValue("");
  cell_finished(a.result, cell);
  if (OutputContainer.items.length > 1) {
    OutputContainer.body.dom.scrollTop =
      OutputContainer.body.dom.scrollHeight -
      cell.body.dom.scrollHeight - 20;
  }
}

var repl_input = new Ext.FormPanel({
  title: 'YT Input',
      url: 'push',
      flex: 0.2,
      layout: 'fit',
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
	  }
	},
	  },],
      });

var NorthButton = new Ext.Button(
				 {text : 'North',
				     pageX : 205,
				     pageY : 10
//        handler: function(b, e) { window.open("session.py", "_top"); }
				     }
				 );

var EastButton = new Ext.Button(
				{text:'East',
				    pageX : 410,
				    pageY : 205
				    }
				);
var SouthButton = new Ext.Button(
				 {text:'South',
				     pageX : 205,
				     pageY : 410
				     }
				 );
var WestButton = new Ext.Button(
				{text:'West',
				    pageX : 10,
				    pageY : 205
				    }
				);

var OutputContainer = new Ext.Panel({
  title: 'YT Output',
      id: 'output_container',
      autoScroll: true,
      flex: 0.8,
      items: []
      });

var PlotPanel = new Ext.Panel(
			      {
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
      title: 'Objects',
      layout: 'anchor',
      region:'west',
      split: true,
      anchor: '100% -35',
      minSize: 150,
      autoScroll: true,
      rootVisible: false,
      root:new Ext.tree.TreeNode({
	expanded:true
	    ,leaf:false
	    ,text:''
	    })
      });

var ButtonGroupPanel = new Ext.Panel({
  layout: 'anchor',
      ButtonAlign: 'center',
      collapsible: false,
      renderTo: document.body,
      tbar: [{
      xtype: 'buttongroup',
	  columns: 3,
	  items: [{
	  text: 'Download',
	      layout:'anchor',
	      anchor: '100% 33%',
	      handler: function(b, e) { window.open("session.py", "_top"); }
	  },{
	  text: 'Save',
	      layout:'anchor',
	      anchor: '100% 67%',
	      },{
	  text: 'Pastebin',
                layout:'anchor',
                anchor: '100% 100%',
		}]
	  }]
      });

var status_panel;
var logging_store = new Ext.data.Store({
  fields: [{name:'record'}],
      reader: new Ext.data.ArrayReader({}, [{name: 'record'}]),
      });

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
		  {
		    // lazily created panel (xtype:'panel' is default)
		  xtype: 'grid',
		      store: logging_store,
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
		      title: 'BETA Sequences',
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
		      },
		  // in this instance the TabPanel is not wrapped by another panel
		  // since no title is needed, this Panel is added directly
		  // as a Container
		  {
		  xtype: 'tabpanel',
		      region: 'center', // a center region is ALWAYS required for border layout
		      id: 'center-panel',
		      deferredRender: false,
		      activeTab: 0,     // first tab initially active
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
			    }, PlotPanel]
		      //                }, {
		      //                   title: 'Plot Window 1',
		      //                  iconCls: 'graph',
		      //                 layout:'anchor',
		      //                anchor:'100% 100%',
		      //               closable: true,
		      //              autoScroll: true,
		      //             items: [ PlotPanel ]
		      //                
		      //       }]
		      }]
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
					      {record: '4d3d3d3 engaged' });
    logging_store.add(record, number_log_records++);
    repl_input.get("input_line").focus();
  });
