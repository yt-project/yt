/**********************************************************************
Debug helper for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.controller.Debugger', {
    extend: 'Ext.app.Controller',
    stores: ['WidgetTypes', 'WidgetInstances', 'Requests'],

    getWidget: function(widgetId) {
        this.getWidgetInstancesStore().find(
            {'widgetid': widgetId});
        var widgetInfo = this.getWidgetInstancesStore().getAt(resultId).data;
        return widgetInfo.widget;
    },

    getAllWidgetsByType: function(typeName) {
        var arr = []
        this.getWidgetInstancesStore().each(function(r) {
            if (r.data['widgettype'] == typeName) {
                arr.push(r.data);
            }
        });
        return arr;
    },

});

