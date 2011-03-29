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
