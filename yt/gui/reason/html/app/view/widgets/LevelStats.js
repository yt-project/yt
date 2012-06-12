/**********************************************************************
Level display stats

Author: Cameron Hummels <chummels@gmail.com>
Affiliation: Columbia
Author: Jeffrey S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
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
        minimum: 0.0,
        grid: true,
        title: 'Relative Portion of Mesh',
        roundToDecimal: function(v) { console.log(v); return(v); },
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
