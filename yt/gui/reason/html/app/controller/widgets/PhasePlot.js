/**********************************************************************
The Plot Window Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/


Ext.define("Reason.controller.widgets.PhasePlot", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.PhasePlotCreator',
               'Reason.view.widgets.PhasePlot'],
    templates: {
        createPhasePlot: 'widget_store.create_phase({objname}, ' + 
                         '"{xField}", "{yField}", "{zField}", {weight})',
        setupPlot:       'widget_store["{widget.varname}"]._setup_plot()',
    },

    executionTriggers: [
        ['#imagepanel', 'afterrender', 'setupPlot'],
    ],

    viewRefs: [
        { ref: 'yTicks', selector: '#y_ticks'},
        { ref: 'xTicks', selector: '#x_ticks'},
        { ref: 'colorTicks', selector: '#cb_ticks'},
        { ref: 'colorbar', selector: '#colorbar'},
        { ref: 'image', selector: '#imagepanel'},
        { ref: 'metadataString', selector: '#metadataString'},
    ],
    
    applyPayload: function(payload) {
        this.getImage().getEl().dom.src = 
            "data:image/png;base64," + payload['image_data'];
        this.getColorbar().getEl().dom.src=
            "data:image/png;base64," + payload['cbar']['cmap_image'];
        var YTicks = this.getYTicks();
        YTicks.removeAll();
        Ext.each(payload['yax']['ticks'], function(tick, index) {
            YTicks.add({xtype:'panel',
                       width: 20, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: right;' +
                              'padding-right: 5px;',
                       html: ' ' + tick[2] + ' ',
                       x:0, y: tick[0]-6});
            YTicks.add({xtype:'panel',
                       width: 20, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:20, y: tick[0]});
        });
        YTicks.doLayout();
        var XTicks = this.getXTicks();
        XTicks.removeAll();
        Ext.each(payload['xax']['ticks'], function(tick, index) {
            XTicks.add({xtype:'panel',
                       width: 1, height:20,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:(400 - tick[0]) + 10, y: 0});
            XTicks.add({xtype:'panel',
                       width: 20, height:20,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px; text-align: center;',
                       html: ' ' + tick[2] + ' ',
                       x: (400 - tick[0]), y: 20});
        });
        XTicks.doLayout();
        var colorTicks = this.getColorTicks();
        colorTicks.removeAll();
        Ext.each(payload['cbar']['ticks'], function(tick, index) {
            colorTicks.add({xtype:'panel',
                       width: 10, height:1,
                       style: 'background-color: #000000;',
                       html:'&nbsp;',
                       x:0, y: tick[0]});
            colorTicks.add({xtype:'panel',
                       width: 30, height:15,
                       border: false,
                       style: 'font-family: "Inconsolata", monospace;' +
                              'font-size: 12px;',
                       html: ' ' + tick[2] + ' ',
                       x:18, y: tick[0]-6});
        });
        colorTicks.doLayout();
    },

    createView: function() {
        var wd = this.payload['data'];
        this.dataView = Ext.widget("phaseplotwindow",{
            varname : this.payload['varname'],
            title: wd['title'],
        });
        this.createMyRefs(this.dataView.id);
        this.applyExecuteHandlers(this.dataView);
        this.getMetadataString().update(this.payload['data']['metadata_string']);
        this.getMetadataString().mds = this.payload['data']['metadata_string'];
        return this.dataView;
    },

    statics: {
        widgetName: 'phase_plot',
        supportsDataObjects: true,
        supportsParameterFiles: false,
        displayName: 'Phase Plot',
        preCreation: function(obj) {
            var widget = Ext.create(this.getName());
            function makePlot(b, e) {
                var conf = {
                    objname: obj.varname,
                    xField: win.query("#x_field")[0].getValue(),
                    yField: win.query("#y_field")[0].getValue(),
                    zField: win.query("#z_field")[0].getValue(),
                };
                var weight = win.query("#weight")[0].getValue();
                if (weight == 'Mass') {
                    weightField = '"CellMassMsun"';
                } else if (weight == 'Volume') {
                    weightField = '"CellVolume"';
                } else if (weight == 'Just Sum') {
                    weightField = 'None';
                }
                conf['weight'] = weightField;
                var cmd = widget.templateManager.applyObject(
                    conf, 'createPhasePlot');
                reason.server.execute(cmd);
                win.destroy();
            }
            win = Ext.widget("phaseplotcreator", {obj:obj});
            var xFieldObj = win.query("#x_field")[0];
            var yFieldObj = win.query("#y_field")[0];
            var zFieldObj = win.query("#z_field")[0];
            this.fieldStore = Ext.create("Reason.store.Fields");
            this.fieldStore.loadData(obj.field_list);
            xFieldObj.bindStore(this.fieldStore);
            yFieldObj.bindStore(this.fieldStore);
            zFieldObj.bindStore(this.fieldStore);
            xFieldObj.setValue("Density");
            yFieldObj.setValue("Temperature");
            zFieldObj.setValue("CellMassMsun");
            win.query("#create")[0].on('click', makePlot);
            win.query("#cancel")[0].on('click', function(){win.destroy();});
            win.show();
        },
    },
});
