/**********************************************************************
Level display stats

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.LevelStats", {
    extend: 'Ext.chart.Chart',
    title: 'Level Stats',
    alias: 'widget.levelstats',
    layout: 'absolute',
    flex: 1.0,
    itemId: 'levelStats',
    animate: true,
    legend: {
       position: 'right',
    },
    axes: [ {
        type: 'Numeric',
        position: 'bottom',
        fields: ['grid_rel', 'cell_rel'],
        minimum: 0,
        grid: true,
        title: 'Relative Portion of Mesh',
    }, {
        type: 'Category',
        position: 'left',
        fields: ['level'],
        title: 'Level',
    } ],
    series: [{
       type: 'bar',
       axis: 'bottom',
       highlight: true,
       xField: 'level',
       yField: ['grid_rel', 'cell_rel'],
       title: ['Level by Grid Count', 'Level by Cell Count'],
       tips: {
         trackMouse: true,
         width: 140,
         height: 64,
         renderer: function(storeItem, item) {
           this.setTitle('Level ' + storeItem.get('level') + ': ' +
                          storeItem.get('grid_count') + ' grids and ' +
                          storeItem.get('cell_count') + ' cells.');
         },
       },
       animate: true,
    }],
});
