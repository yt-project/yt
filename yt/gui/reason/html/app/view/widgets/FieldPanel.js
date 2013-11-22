/**********************************************************************
Field Info Display

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.FieldPanel", {
    extend: 'Ext.panel.Panel',
    alias: 'widget.fieldpanel',
    iconCls: 'graph',
    autoScroll: true,
    layout: {
        type: 'vbox',
        pack: 'start',
        align: 'stretch',
    },
    width: '100%',
    height: '100%',

    items: [ 
        {
          xtype: 'combo',
          text: 'Field',
          itemId: 'fieldSelector',
          height: 30,
          width: 700,
          queryMode: 'local',
          editable: false,
          triggerAction: 'all',
          validateOnBlur: false,
        }, {
          xtype: 'panel',
          title: 'Field Source',
          itemId: 'fieldSourcePanel',
          flex: 3.0,
          width: 700,
          autoScroll: true,
          bodyCls: 'pfdisplay',
        }, {
          xtype: 'panel',
          title: 'Field Parameters',
          height: 200,
          width: 700,
          autoScroll: true,
          bodyCls: 'pfdisplay',
        }
    ],
});

