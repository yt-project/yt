/**********************************************************************
The Pannable Map Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.view.widgets.PannableMapView", {
    extend: 'Ext.panel.Panel',
    title: 'This should not be visible.',
    alias: 'widget.mapview',
    iconCls: 'graph',
    autoScroll: true,
    layout: 'absolute',
    width: '100%',
    height: '100%',
    closable: true,
    
    items: [ 
        {
          xtype:'panel',
          layout: 'absolute',
          x: 0,
          y: 0,
          width: 512,
          height: 512,
          itemId: 'mapbox',
        }  
    ],
});
