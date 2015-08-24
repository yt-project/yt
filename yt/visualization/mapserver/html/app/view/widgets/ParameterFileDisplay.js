/**********************************************************************
The Parameter File Display Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.ParameterFileDisplay", {
    extend: 'Ext.tab.Panel',
    title: 'This should not be visible.',
    alias: 'widget.pfdisplay',
    iconCls: 'graph',
    autoScroll: true,
    layout: 'absolute',
    width: '100%',
    height: '100%',
    closable: true,
    activeTab: 0,

    items: [
        { 
          xtype: 'panel',
          width: 600,
          height: 500,
          title: 'Dataset Information',
          layout: {
            type: 'absolute',
          },
          items: [
            {
              xtype: 'panel',
              itemId: 'pfParams',
              bodyCls: 'pfdisplay',
              width: 600,
              height: 300,
              x: 0,
              y: 0,
              items: [],
            }, {
              xtype: 'panel',
              itemId: 'widgetpanel',
              width: 600,
              height: 300,
              x: 0,
              y: 300,
              items: [],
            }
          ],
        }, {
          xtype: 'panel',
          itemId: 'fieldPanel',
          title: 'Field Info',
          layout: {
            type: 'hbox',
            pack: 'start',
            align: 'stretch',
          },
          items: []
        }, {
          xtype: 'panel',
          itemId: 'statsPanel',
          title: 'Mesh Statistics',
          layout: {
            type: 'hbox',
            pack: 'start',
            align: 'stretch',
          },
          items: []
        }
    ],
});

