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
               "Reason.controller.widgets.ProgressBar",
    ],
    stores: ['WidgetTypes', 'WidgetInstances'],
    views: ['WidgetTypesGrid', 'WidgetInstancesGrid'],

    init: function() {
        Ext.iterate(Reason.controller.widgets, function(i, w, ws) {
            Ext.require(w.getName());
            this.registerWidget(w);
        }, this);
        this.application.addListener({
            createwidget: {fn: this.createWidget, scope: this},
            showwidgets: {fn: this.showWidgetMenu, scope: this},
            payloadwidget: {fn: this.newWidgetCreated, scope: this},
            payloadwidget_payload: {fn: this.sendPayload, scope: this},
            enabledebug: {fn: this.enableDebug, scope: this},
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

    newWidgetCreated: function(payload) {
        /* We have the following fields:
                type             ('widget')
                widget_type
                varname
                data             (has subfields)

           We now obtain our class, create that with the factory, and we add
           the resultant class to our store.
        */
        var resultId = this.getWidgetTypesStore().find(
            'widgetname', payload['widget_type']);
        if (resultId == -1) {
            Ext.Error.raise('Did not recognize widget type "' +
                            payload['widget_type'] + '".');
        }
        var widgetInfo = this.getWidgetTypesStore().getAt(resultId).data;
        /* The widget adds its view to the viewport. */
        var newWidget = Ext.create(widgetInfo['widgetclass'].getName(),
                            {payload: payload});
        console.log("Adding widget payload with varname " + payload['varname']);
        this.getWidgetInstancesStore().add({
            widgetid: payload['varname'],
            widgettype: widgetInfo.widgetname,
            widget: newWidget
        });
        newWidget.createView();
    },

    sendPayload: function(payload) {
        var resultId = this.getWidgetInstancesStore().find(
            'widgetid', payload['widget_id']);
        examine = payload;
        if (resultId == -1) {
            Ext.Error.raise('Could not find widget "' +
                            payload['widget_id'] + '".');
        }
        var widgetInfo = this.getWidgetInstancesStore().getAt(resultId).data;
        widgetInfo['widget'].applyPayload(payload);
    },

    enableDebug: function() {
        if(this.instanceView) {return;}
        this.instanceView = Ext.widget('widgetinstancesgrid');
        this.typeView = Ext.widget('widgettypesgrid');
        Ext.ComponentQuery.query("viewport > #center-panel")[0].add(
            this.instanceView);
        Ext.ComponentQuery.query("viewport > #center-panel")[0].add(
            this.typeView);
    },

});
