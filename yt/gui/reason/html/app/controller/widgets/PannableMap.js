/**********************************************************************
The Pannable Map Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.controller.widgets.PannableMap", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: ['Reason.view.widgets.PannableMapView'],

    templates: {
        title: "Map View for {widget.varname}",
    },

    widgetTriggers: [
        ['#mapbox', 'afterlayout', 'setupLeaflet'],
    ],

    executionTriggers: [
    ],

    viewRefs: [
        { ref: 'mapBox', selector: '#mapbox' },
    ],

    applyPayload: function(payload) {
        return;
    },

    setupLeaflet: function() {
        var toRender = this.getMapBox().getEl().dom.childNodes[0].id;
        this.leafletMap = new L.Map(toRender,
            {
                center: new L.LatLng(0, 0),
                zoom: 1,
            });
        var ytMapURL = this.payload.data['prefix'] + '/map/{z}/{x}/{y}.png';
        this.tileLayer = new L.TileLayer(ytMapURL, {maxZoom: 18});
        this.leafletMap.addLayer(this.tileLayer)
    },

    createView: function() {
        this.dataView = Ext.widget("mapview", {
             title: 'Pannable Map of ' + this.payload['field'],
        });
        this.createMyRefs(this.dataView.id);
        this.applyExecuteHandlers(this.dataView);
        return this.dataView;
    },

    statics: {
        widgetName: "pannable_map",
        displayName: "Pannable Map",
        supportsDataObjects: false,
        supportsParameterFiles: false,
        preCreation: function(obj) { },
    },
});
