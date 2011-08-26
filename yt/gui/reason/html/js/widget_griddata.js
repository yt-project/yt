/**********************************************************************
The Plot Window Widget

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

// shim layer with setTimeout fallback
var WidgetGridDataViewer = function(python_varname, widget_data) {
    this.id = python_varname;
    this.widget_data = widget_data;
    store = new Ext.data.ArrayStore({
        fields: [
           {name: 'grid_id', type:'int'},
           {name: 'level', type:'int'},
           {name: 'left_edge_x', type: 'float'},
           {name: 'left_edge_y', type: 'float'},
           {name: 'left_edge_z', type: 'float'},
           {name: 'right_edge_x', type: 'float'},
           {name: 'right_edge_y', type: 'float'},
           {name: 'right_edge_z', type: 'float'},
           {name: 'dim_x', type: 'int'},
           {name: 'dim_y', type: 'int'},
           {name: 'dim_z', type: 'int'},
           {name: 'cells', type: 'int'},
        ]
    });
    store.loadData(widget_data['gridvals']);
    examine = widget_data;

    viewport.get("center-panel").add(
        {
            xtype: 'panel',
            id: "gg_" + python_varname,
            title: "Grid Data Viewer",
            iconCls: 'graph',
            autoScroll: true,
            layout:'vbox',
            layoutConfig: {align: 'stretch', pack: 'start'},
            closable: true,
            items: [ {
                       xtype: 'grid',
                       store: store,
                       columns: [
                            {
                                id: 'grid_id',
                                header: 'Grid ID',
                                width: 100,
                                dataIndex: 'grid_id',
                                sortable: true,
                            }, {
                                id: 'left_edge_x',
                                header: 'Left Edge x',
                                width: 100,
                                dataIndex: 'left_edge_x',
                                sortable: true,
                            }, {
                                id: 'left_edge_y',
                                header: 'Left Edge y',
                                width: 100,
                                dataIndex: 'left_edge_y',
                                sortable: true,
                            }, {
                                id: 'left_edge_z',
                                header: 'Left Edge z',
                                width: 100,
                                dataIndex: 'left_edge_z',
                                sortable: true,
                            }, {
                                id: 'right_edge_x',
                                header: 'Right Edge x',
                                width: 100,
                                dataIndex: 'right_edge_x',
                                sortable: true,
                            }, {
                                id: 'right_edge_y',
                                header: 'Right Edge y',
                                width: 100,
                                dataIndex: 'right_edge_y',
                                sortable: true,
                            }, {
                                id: 'right_edge_z',
                                header: 'Right Edge z',
                                width: 100,
                                dataIndex: 'right_edge_z',
                                sortable: true,
                            }, {
                                id: 'dim_x',
                                header: 'DimX',
                                width: 100,
                                dataIndex: 'dim_x',
                                sortable: true,
                            }, {
                                id: 'dim_y',
                                header: 'DimY',
                                width: 100,
                                dataIndex: 'dim_y',
                                sortable: true,
                            }, {
                                id: 'dim_z',
                                header: 'DimZ',
                                width: 100,
                                dataIndex: 'dim_z',
                                sortable: true,
                            }, {
                                id: 'cells',
                                header: 'Cells',
                                width: 100,
                                dataIndex: 'cells',
                                sortable: true,
                            },
                       ],
                      flex: 1,
                      }
                   ],

        }
    );

    viewport.get("center-panel").activate("gg_" + this.id);
    viewport.doLayout();
    this.panel = viewport.get("center-panel").get("gg_" + python_varname);
    this.panel.doLayout();

    this.accept_results = function(payload) { }
}

widget_types['grid_data'] = WidgetGridDataViewer;
