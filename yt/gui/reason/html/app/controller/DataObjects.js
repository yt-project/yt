/**********************************************************************
Data object controller for Reason

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
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
           payloaddataobjects : {fn: this.refreshDataObjects, scope: this},
        });
        this.control({
            "#dataobjects": { itemcontextmenu:
                function(view, rec, node, index, e) {
                    e.preventDefault();
                    this.application.fireEvent("showwidgets", rec, e);
            }
        }});
        this.callParent(arguments);
    },

    refreshDataObjects: function(payload) {
        /*console.log("Refreshing data objects");*/
        var objs = payload['objs'];
        var root = this.getDataObjectsStore().getRootNode();
        root.removeAll();
        var pf;
        Ext.each(objs, function(o, i, os) {
            /*console.log("Appending " + o['name']);*/
            pf = root.appendChild({
                name: o.name,
                type: o.type,
                filename: o.filename,
                field_list: o.field_list,
                varname: o.varname,
                leaf: false,
                expanded: true,
                iconCls: 'pf_icon'
            });
            Ext.each(o['children'], function(c, ci, cs) {
                /*console.log("    Appending " + c['name']);*/
                pf.appendChild({name: c.name,
                                type: c.type,
                                filename: o.filename,
                                field_list: o.field_list,
                                varname: c.varname,
                                leaf: true,
                                iconcls: 'data_obj'});

            });
        });
    },
});

