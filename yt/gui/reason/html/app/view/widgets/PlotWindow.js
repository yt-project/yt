/**********************************************************************
The Plot Window Widget View

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

Ext.define("Reason.view.widgets.PlotWindow", {
    extend: 'Ext.panel.Panel',
    title: 'This should not be visible.',
    alias: 'widget.plotwindow',
    iconCls: 'graph',
    autoScroll: true,
    layout:'absolute',
    width: '100%',
    height: '100%',
    closable: true,

    items: [ 
        { xtype:'panel',
          id: 'image_panel',
          autoEl: {
                tag: 'img',
                id: "mainImage",
                width: 400,
                height: 400,
                draggable: false,
            },
            x: 100,
            y: 10,
            width: 400,
            height: 400,
            draggable: false,
        }, {
            xtype:'panel',
            id: 'colorbar',
            autoEl: {
                tag: 'img',
                id: 'cb',
                src: '',
                width: 28,
                height: 398,
                style: 'border: 1px solid #000000;',
            },
            x: 510,
            y: 10,
            width: 30,
            height: 400,
        }, {
            xtype: 'panel',
            id: 'ticks',
            layout: 'absolute',
            y: 0,
            x: 540,
            width: 100,
            height: 420,
            items : [],
            border: false,
        }, {   xtype: 'multislider',
            id: 'slider',
            minValue: 0,
            maxValue: 100,
            increment: 0.1,
            x: 100, y: 410,
            width: 400,
        },{
            xtype: 'combo',
            text: 'Field',
            x: 100,
            y: 435,
            width: 400,
            store: [],
            value: '',
            editable: false,
            triggerAction: 'all',
            validateOnBlur: false,
        }, {
        /* the single buttons for 10% pan*/
            xtype:'button',
            iconCls: 'singleuparrow',
            id: 'singleuparrow',
            //text: 'North',
            x: 40,
            y: 10,
        }, {
            xtype:'button',
            iconCls: 'singlerightarrow',
            id: 'singlerightarrow',
            //text:'East',
            x : 60,
            y : 30,
        }, {
            xtype:'button',
            iconCls: 'singledownarrow',
            id: 'singledownarrow',
            //text: 'South',
            x: 40,
            y: 50,
        }, {
            xtype: 'button',
            iconCls: 'singleleftarrow',
            id: 'singleleftarrow',
            //text: 'West',
            x: 20,
            y: 30,
        }, 
        /* the double buttons for 50% pan*/
        {
            xtype:'button',
            iconCls: 'doubleuparrow',
            id:'doubleuparrow',
            //text: 'North',
            x: 40,
            y: 80,
        }, {
            xtype:'button',
            iconCls: 'doublerightarrow',
            id:'doublerightarrow',
            //text:'East',
            x : 60,
            y : 100,
        }, {
            xtype:'button',
            iconCls: 'doubledownarrow',
            //text: 'South',
            id: 'doubledownarrow',
            x: 40,
            y: 120,
        }, {
            xtype: 'button',
            iconCls: 'doubleleftarrow',
            id: 'doubleleftarrow',
            //text: 'West',
            x: 20,
            y: 100,
        },
        /* Now the zoom buttons */
        {
            xtype: 'button',
            text: 'Zoom In 10x',
            id: "zoom10x",
            x: 10,
            y: 160,
            width: 80,
        },{
            xtype: 'button',
            text: 'Zoom In 2x',
            id: "zoom2x",
            x: 10,
            y: 185,
            width: 80,
        },{
            xtype: 'button',
            text: 'Zoom Out 2x',
            id:'zoomout2x',
            x: 10,
            y: 210,
            width: 80,
        },{
            xtype: 'button',
            text: 'Zoom Out 10x',
            id:'zoomout10x',
            x: 10,
            y: 235,
            width: 80,
        },{
            xtype: 'button',
            text: 'Upload Image',
            x: 10,
            y: 285,
            width: 80,
            tooltip: "Upload the current image to " +
                     "<a href='http://imgur.com'>imgur.com</a>",
        },{
            xtype: 'button',
            text: 'Pannable Map',
            x: 10,
            y: 335,
            width: 80,
            tooltip: "Open a pannable map in a new tab",
        },{
            xtype: 'panel',
            layout: 'vbox',
            id: 'rhs_panel',
            width: 250,
            height: 460,
            x: 690, y: 10,
            layoutConfig: {
                align: 'stretch',
                pack: 'start',
            },
            items: [
                {
                  xtype: 'panel',
                  title: 'Plot MetaData',
                  id: 'metadata',
                  style: {fontFamily: '"Inconsolata", monospace'},
                  html: 'Welcome to the Plot Window.',
                  height: 200,
                }, {
                  xtype: 'tabpanel',
                  id: 'editor_panel',
                  flex: 1,
                  activeTab: 0,
                  items: [
                {
                  xtype: 'panel',
                  title: 'Plot Editor',
                  id: 'plot_edit',
                  style: {fontFamily: '"Inconsolata", monospace'},
                  layout: 'absolute',
                  flex: 1,
                  items : [
                     {
                       x: 10,
                       y: 20,
                       width: 70,
                       xtype: 'label',
                       text: 'Display',
                     },
                     {
                       x: 80,
                       y: 20,
                       width : 80,
                       xtype: 'combo',
                       editable: false,
                       triggerAction: 'all',
                       validateOnBlur: false,
                       store: ['log10', 'linear'],
                       value: 'linear',
                     },
                     {
                       x: 10,
                       y: 60,
                       width: 70,
                       xtype: 'label',
                       text: 'Colormap',
                     },
                     {
                       x: 80,
                       y: 60,
                       width : 140,
                       xtype: 'combo',
                       editable: false,
                       triggerAction: 'all',
                       validateOnBlur: false,
                       store: ['algae', 'RdBu', 'gist_stern',  
                               'hot', 'jet', 'kamae', 
                                'B-W LINEAR', 'BLUE',
                                'GRN-RED-BLU-WHT', 'RED TEMPERATURE',
                                'BLUE', 'STD GAMMA-II', 'PRISM',
                                'RED-PURPLE', 'GREEN', 'GRN',
                                'GREEN-PINK', 'BLUE-RED', '16 LEVEL',
                                'RAINBOW', 'STEPS', 'STERN SPECIAL',
                                'Haze', 'Blue - Pastel - Red',
                                'Pastels', 'Hue Sat Lightness 1',
                                'Hue Sat Lightness 2', 'Hue Sat Value 1',
                                'Hue Sat Value 2', 'Purple-Red + Stripes',
                                'Beach', 'Mac Style', 'Eos A', 'Eos B',
                                'Hardcandy', 'Nature', 'Ocean', 'Peppermint',
                                'Plasma', 'Blue-Red', 'Rainbow', 'Blue Waves',
                                'Volcano', 'Waves', 'Rainbow18',
                                'Rainbow + white', 'Rainbow + black'],
                       value: 'algae',
                     }
                  ]
                }, {
                  xtype: 'panel',
                  title: 'Contours',
                  id: 'contour_edit',
                  style: {fontFamily: '"Inconsolata", monospace'},
                  layout: 'absolute',
                  flex: 1,
                  items : [
                     {
                       x: 10,
                       y: 20,
                       width: 70,
                       xtype: 'label',
                       text: 'Field',
                     }, {
                       x: 80,
                       y: 20,
                       width : 160,
                       xtype: 'combo',
                       editable: false,
                       id: 'field',
                       triggerAction: 'all',
                       validateOnBlur: false,
                       value: '',
                       store: [],
                     }, {
                       x: 10,
                       y: 60,
                       width: 70,
                       xtype: 'label',
                       text: 'Levels',
                     }, {
                       x: 80,
                       y: 60,
                       width : 160,
                       xtype: 'slider',
                       id: 'ncont',
                       minValue: 0,
                       maxValue: 10,
                       value: 5,
                       increment: 1,
                       /*plugins: [ {ptype: 'slidertip'} ],*/
                     }, {
                       x: 10,
                       y: 100,
                       width: 70,
                       xtype: 'label',
                       text: 'Logspaced',
                     }, {
                       x: 80,
                       y: 100,
                       width : 160,
                       xtype: 'checkbox',
                       id: 'logit',
                       checked: true,
                     }, {
                       x: 10,
                       y: 180,
                       width: 80,
                       xtype: 'button',
                       text: 'Apply',
                     }
                  ]
                }, {
                  xtype: 'panel',
                  title: 'Velocity Vectors',
                  id: 'vector_edit',
                  style: {fontFamily: '"Inconsolata", monospace'},
                  layout: 'absolute',
                  flex: 1,
                  items : [
                     {
                       x: 10,
                       y: 60,
                       width: 70,
                       xtype: 'label',
                       text: 'Skip Factor',
                     }, {
                       x: 80,
                       y: 60,
                       width : 160,
                       xtype: 'slider',
                       id: 'skip',
                       minValue: 1,
                       maxValue: 64,
                       value: 32,
                       increment: 1,
                       /*plugins: [ {ptype: 'slidertip'} ],*/
                     }, {
                       x: 10,
                       y: 180,
                       width: 80,
                       xtype: 'button',
                       text: 'Apply',
                     }
                  ]
                }
                ] } /* tabpanel items and entry */
                ]
        }
    ],
});

