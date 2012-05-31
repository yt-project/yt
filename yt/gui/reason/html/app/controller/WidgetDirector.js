/**********************************************************************
Widget controller class

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

Ext.define('Reason.controller.WidgetDirector', {
    extend: 'Ext.app.Controller',
    requires: ["Reason.controller.widgets.SampleWidget",
               "Reason.controller.widgets.PlotWindow",
    ],
    stores: ['WidgetTypes', 'WidgetInstances'],

    init: function() {
        Ext.iterate(Reason.controller.widgets, function(i, w, ws) {
            Ext.require(w.getName());
            this.registerWidget(w);
        }, this);
        this.application.addListener({
            createwidget: {fn: this.createWidget, scope: this},
            showwidgets: {fn: this.showWidgetMenu, scope: this},
        });
        this.callParent(arguments);
    },

    registerWidget: function(w) {
        if (w.widgetName == null) {return;}
        console.log("Registering " + w.widgetName);
        this.getWidgetTypesStore().add({
                   widgetname: w.widgetName,
                   widgetclass: w,
                   displayname: w.displayName,
                   pfs: w.supportsParameterFiles,
                   objs: w.supportsDataObjects,
        });
    },

    createWidget: function(b, e) {
        var w = b.widget;
        console.log("Asked to create " + b.widget.widgetName);
        var store = this.getWidgetInstancesStore();
        b.widget.preCreation(b.dataObj);
    },

    showWidgetMenu: function(treerecord, e) {
        var contextMenu = Ext.create('Ext.menu.Menu', {plain: true,});
        var data = treerecord.data;
        var w;
        this.getWidgetTypesStore().each(function(record, idx) {
            w = record.data;
            examine = w;
            if (((data.type == 'parameter_file') && (w.pfs  == false)) 
             || ((data.type != 'parameter_file') && (w.objs == false))) {
              return;
            }
            contextMenu.add({xtype:'menuitem',
                             text: w.displayname,
                             listeners: {
                                click: {
                                    fn : this.createWidget,
                                    scope: this
                                },
                             },
                             widget: w.widgetclass,
                             dataObj: data
            });
        }, this);
        contextMenu.showAt(e.getXY());
    },

});
