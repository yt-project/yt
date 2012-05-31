/**********************************************************************
Data object controller for Reason

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

Ext.define('Reason.controller.DataObjects', {
    extend: 'Ext.app.Controller',
    stores: [ 'DataObjects' ],
    models: ['DataObject'],
    view: ['DataObjectTree'],
    refs: [
        { ref: 'dataObjectTree',
          selector: 'dataobjecttree#dataobjects',
          xtype: 'dataobjecttree',
          itemId: 'dataobjects',
          autoCreate: true,
        },
    ],


    init: function() {
        this.application.addListener({
           newdataobjects : {fn: this.refreshDataObjects, scope: this},
        });
        this.callParent(arguments);
    },

    refreshDataObjects: function(objs) {
        console.log("Refreshing data objects");
        var view = this.getDataObjectTree();
        var store = this.getDataObjectsStore();
        store.removeAll();
        store.setRootNode({text: '', leaf: false, expanded: true});
        var root = store.getRootNode();
        var pf;
        Ext.each(objs, function(o, i, os) {
            console.log("Appending " + o['name']);
            pf = root.appendChild({
                text: o['name'],
                objdata: o,
                leaf: false,
                expanded: true,
                iconCls: 'pf_icon'
            });
            Ext.each(o['children'], function(c, ci, cs) {
                console.log("    Appending " + c['name']);
                pf.appendChild({text: c['name'],
                                objdata: c,
                                leaf: true,
                                iconcls: 'data_obj'});

            });
        });
    },
});

