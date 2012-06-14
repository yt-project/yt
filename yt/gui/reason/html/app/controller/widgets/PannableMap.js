/**********************************************************************
The Pannable Map Widget

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
        examine = toRender;
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
