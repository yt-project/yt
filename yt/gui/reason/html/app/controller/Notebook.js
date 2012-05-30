/**********************************************************************
Notebook controller for Reason

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

Ext.define('Reason.controller.Notebook', {
    extend: 'Ext.app.Controller',
    stores: [ 'CellValues' ],
    views: ['Notebook'],
    refs: [
        { ref: 'inputLine',
          selector: '#input_line'
        }
    ],

    init: function() {
        this.application.addListener({
            newcell: {fn: this.addCell, scope: this},
            executecell: {fn: this.executeCell, scope: this},
        })
        this.control({
            '#executecellbutton': {
                click: function(f, e) {
                    this.executeCell(this.getInputLine().getValue());
                }
            },
            '#inputline': {
                specialkey: function(ed, field, e, opts){
                    if (e.getKey() == e.ENTER) {
                        this.executeCell(field.getValue());
                    }
                },
            },
        });
        this.callParent(arguments);
    },

    addCell: function(cell) {
        examine = this;
        this.getCellValuesStore().add({
            input: cell['input'],
            output: cell['output'],
            raw_input: cell['raw_input'],
            executiontime: cell['executiontime'],
        });
    },
    executeCell: function(line) {
        console.log("Asked to execute " + line);
    }
});

