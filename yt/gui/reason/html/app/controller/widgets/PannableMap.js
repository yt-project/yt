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

var WidgetPannableMap = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "pm_" + this.id,
            title: "Pannable Map",
            iconCls: 'graph',
            autoScroll: true,
            layout:'absolute',
            closable: true,
            items: [ 
                {
                    xtype:'box',
                    autoEl: {
                        tag: 'div',
                        id: "map_" + this.id,
                        width: 512,
                        height: 512,
                    },
                    x: 10,
                    y: 10,
                    width: 512,
                    height: 512,
                    listeners: {afterrender:
                        function() {
                          var map = new L.Map('map_' + python_varname, {
                                  center: new L.LatLng(0.0, 0.0),
                                  zoom: 0,
                                  });
                          var cloudmadeUrl = widget_data['prefix'] + '/map/{z}/{x}/{y}.png';
                          cloudmade = new L.TileLayer(cloudmadeUrl, {maxZoom: 18});
                          map.addLayer(cloudmade);
                    }},
                }  
            ]
        }
    );

    viewport.get("center-panel").activate("pm_" + this.id);
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("pm_" + this.id);
    this.panel.doLayout();
    examine = this.panel;

    this.accept_results = function(payload) { }
}

widget_types['pannable_map'] = WidgetPannableMap;
