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
    requires: ["Reason.controller.widgets.SampleWidget"],
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
        if (w.prototype.widgetName == null) {return;}
        console.log("Registering " + w.prototype.widgetName);
        this.getWidgetTypesStore().add({
                   widgetname: w.prototype.widgetName,
                   widgetclass: w,
                   displayname: w.prototype.displayName,
                   pfs: w.prototype.supportsParameterFiles,
                   objs: w.prototype.supportsDataObjects,
        });
    },

    createWidget: function(b, e) {
        var w = b.widget;
        console.log("Asked to create " + b.widget.widgetName);
        b.widget.createWidget(b.dataObj);
    },

    showWidgetMenu: function(treerecord, e) {
        var contextMenu = Ext.create('Ext.menu.Menu', {plain: true,});
        
        var w;
        examine = treerecord;
        this.getWidgetTypesStore().each(function(record, idx) {
            w = record.data;
            contextMenu.add({xtype:'menuitem',
                             text: w.displayname,
                             handler: this.createWidget,
                             widget: w.widgetclass,
                             dataObj: treerecord.data,
            });
        }, this);
        contextMenu.showAt(e.getXY());
    },

});
