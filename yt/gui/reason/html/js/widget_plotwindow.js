var WidgetPlotWindow = function(python_varname, viewport) {
    this.vn = "pw_" + python_varname;
    this.print_python = function(b, e) {
        yt_rpc.ExtDirectREPL.execute({code:'print 1,2,3, "' + python_varname + '"'},
                                     function(f, a)
                                     {alert(a.result['output']);});
    }

    viewport.get("center-panel").add(
                  {
                  xtype: 'panel',
                  id: this.vn,
                  title: this.vn,
                  iconCls: 'graph',
                  autoScroll: true,
                  layout:'absolute',
                  items: [ 
                      {xtype:'button',
                       text: 'North',
                       x: 205,
                       y: 10},
                      {xtype:'button',
                       text:'East',
                       x : 410,
                       y : 205,
                       handler: this.print_python},
                      {xtype:'button',
                       text: 'South',
                       x: 205,
                       y: 410},
                      {xtype: 'button',
                       text: 'West',
                       x: 10,
                       y: 205},
                       ]
                  });

    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("pw_" + python_varname);
    this.panel.doLayout();
}
