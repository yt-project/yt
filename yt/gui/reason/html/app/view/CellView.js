/**********************************************************************
Notebook Cell view for Reason

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

var cellDisplay = new Ext.XTemplate(
    '<b>Input:</b><br/><br/>',
    '{input}',
    '<br/><br/>',
    '<hr>',
    '<b>Output:</b><br/><br/>',
    '<pre>{output}</pre><br/><br/>'
);

Ext.define('Reason.view.CellView', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.notebookcells',
    title: 'Cells',
    store: 'CellValues',
    itemId: 'cells',
    autoscroll: true,
    flex: 0.7,
    columns: [{header:'Execution Time', dataIndex: 'executiontime', flex:1}],
    viewConfig: {
        stripeRows: false,
        disableSelection: true,
        trackOver: false,
    },
    features: [{
        ftype: 'rowbody',
        getAdditionalData: function(data, rowIndex, record, orig) {
            return {
                rowBody: cellDisplay.apply(data),
                rowBodyCls: 'codeview',
                rowBodyColspan: this.view.headerCt.getColumnCount(),
            };
        }
      }, {
        ftype: 'rowwrap'
      }],
});

