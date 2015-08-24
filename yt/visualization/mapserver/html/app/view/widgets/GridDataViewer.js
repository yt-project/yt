/**********************************************************************
The Grid Data Viewer View

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.GridDataViewer", {
    extend: 'Ext.grid.Panel',
    title: 'This should not be visible.',
    alias: 'widget.griddataviewer',
    iconCls: 'graph',
    autoScroll: true,
    width: '100%',
    height: '100%',
    closable: true,
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
});

