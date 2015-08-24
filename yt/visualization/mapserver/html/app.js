/**********************************************************************
The main GUI facility for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.Loader.setConfig({enabled:true});

var reason, rdebug, examine;

Ext.application({
    requires: ['Ext.container.Viewport',
               'Reason.view.LoggingGrid',
               'Reason.view.DataObjectTree',
               'Reason.controller.widgets.SampleWidget',
               ],
    name: 'Reason',

    appFolder: 'reason/app',

    launch: function() {
        reason = this;
        this.numberOfRequests = 0;
        rdebug = this.getController("Debugger");
        this.server = this.getController("ServerCommunication");
        Ext.create('Ext.container.Viewport', {
            layout: 'border',
            items: [
                {
                    xtype: 'logginggrid',
                    autofill: true,
                    region: 'south',
                    id: "status-region",
                    cls: "status-logger",
                    split: true,
                    height: 100,
                    maxSize: 200,
                    collapsible: true,
                    title: 'Status',
                    margins: '0 0 0 0',
               }, {
                    id: 'west-panel',
                    region: 'west',
                    split: true,
                    width: 200,
                    minSize: 175,
                    maxSize: 400,
                    collapsible: true,
                    layout: {
                        type: 'vbox',
                        pack: 'start',
                        align: 'stretch',
                    },
                    items: [
                        { xtype: 'dataobjecttree', 
                          flex: 1.0},
                        { xtype: 'panel',
                          tpl: 'Pending Requests: {0}',
                          data: [0],
                          height: 30,
                          id: 'pending',
                        }
                    ],
              }, {
                xtype: 'tabpanel',
                region: 'center',
                id: 'center-panel',
                activeTab: 0,
                dockedItems: [ {
                    dock: 'top',
                    xtype: 'mainmenu',
                } ],
                items: [
                { xtype: 'panel',
                  title: 'Welcome!',
                  closeable: true,
                  autoLoad: 'help.html',
                },
                { xtype: 'ytnotebook' },
                ]
              }
            ]
        });
        this.pending = Ext.getCmp("pending");
        this.fireEvent("allowinput");
    },
    controllers : [
        'Logging',
        'DataObjects',
        'Notebook',
        'ServerCommunication',
        'WidgetDirector',
        'MenuActions',
        'Debugger',
        'FileOpen',
    ],

});
