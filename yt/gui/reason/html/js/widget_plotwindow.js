var WidgetPlotWindow = function(python_varname) {
    this.id = python_varname;
    this.print_python = function(b, e) {
        yt_rpc.ExtDirectREPL.execute(
            {code:'print "' + python_varname + '"'},
            function(f, a) {alert(a.result['output']);}
        );
    }

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "pw_" + this.id,
            title: this.id,
            iconCls: 'graph',
            autoScroll: true,
            layout:'absolute',
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
                    x: 10,
                    y: 10,
                    width: 400,
                    height: 400,
                }, {
                    xtype:'button',
                    text: 'North',
                    x: 205,
                    y: 10,
                    handler: function(b,e) {
                        cc = python_varname + '.pan_rel((0.0, 0.5))'
                        yt_rpc.ExtDirectREPL.execute(
                        {code:cc}, handle_payload); 
                    }
                }, {
                    xtype:'button',
                    text:'East',
                    x : 410,
                    y : 205,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.5, 0.0))'},
                        handle_payload); 
                    }
                }, {
                    xtype:'button',
                    text: 'South',
                    x: 205,
                    y: 410,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((0.0, -0.5))'},
                        handle_payload); 
                    }
                }, {
                    xtype: 'button',
                    text: 'West',
                    x: 10,
                    y: 205,
                    handler: function(b,e) {
                        yt_rpc.ExtDirectREPL.execute(
                            {code:python_varname + '.pan_rel((-0.5, 0.0))'},
                        handle_payload); 
                    }
                },
            ]
        }
    );

    viewport.get("center-panel").activate("pw_" + this.id);
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("pw_" + python_varname);
    this.panel.doLayout();
    this.image_panel = this.panel.get("image_panel_"+python_varname);

    this.accept_results = function(payload) {
        this.image_panel.el.dom.src = "data:image/png;base64," + payload['image_data'];
    }

    yt_rpc.ExtDirectREPL.execute(
        {code:python_varname + '.zoom(1.0)'},
        handle_payload);
}

widget_types['plot_window'] = WidgetPlotWindow;
