/**********************************************************************
The Plot Window Widget View

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

Ext.define("Reason.view.widgets.Scene", {
    extend: 'Ext.panel.Panel',
    title: 'This should not be visible.',
    alias: 'widget.scene',
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
          x: 10,
          y: 10,
          width: 450,
          height: 450,
          itemId: 'scenepanel',
        }, {
            xtype: 'multislider',
            itemId: 'cameraPathSlider',
            minValue: 0,
            maxValue: 100,
            increment: 1,
            x: 10, y: 460,
            width: 450,
        },{
            xtype: 'panel',
            itemId: 'rhs_panel',
            width: 300,
            height: 450,
            x: 470, y: 10,
            layout: {
                type: 'vbox',
                align:'stretch',
                pack:'start',
            },
            items: [
                {
                  xtype: 'tabpanel',
                  title: 'Scene Editor',
                  itemId: 'editor_panel',
                  flex: 1.0,
                  activeTab: 0,
                  items: [
                {
                  xtype: 'panel',
                  title: 'Widgets',
                  itemId: 'widget_edit',
                  flex: 1,
                  layout: {
                    type: 'vbox',
                    pack: 'start',
                    align: 'stretch',
                  },
                  items : [
                    {
                      xtype: 'gridpanel',
                      itemId: 'widgetlist',
                      flex: 1.0,
                      selType: 'rowmodel',
                      columns: [
                        {
                          itemId: 'widgetEnabled',
                          header: '',
                          dataIndex: 'enabled',
                          xtype: 'checkcolumn',
                          width: 30,
                        }, {
                          itemId: 'name',
                          header: 'Name',
                          dataIndex: 'name',
                          editor: false,
                        }, {
                          itemId: 'type',
                          header: 'Type',
                          dataIndex: 'type',
                          editor: false,
                        }
                      ],
                      plugins: [
                        {ptype: 'cellediting',
                         clicksToEdit: 1},
                      ],
                    },
                  ]
                }, {
                  xtype: 'panel',
                  title: 'Camera Path',
                  itemId: 'camera_path',
                  layout: {
                    type: 'vbox',
                    pack: 'start',
                    align: 'stretch',
                  },
                  flex: 1.0,
                  items : [
                    {
                      xtype: 'gridpanel',
                      itemId: 'keyframeview',
                      flex: 1.0,
                      selType: 'rowmodel',
                      columns: [
                        {
                          id: 'time',
                          header: 'Time',
                          dataIndex: 'time',
                          editor: {
                            xtype: 'textfield',
                          },
                          flex: 1.0,
                          sortable: true,
                        }, {
                          id: 'pos_x',
                          header: 'Pos X',
                          dataIndex: 'pos_x',
                          editor: false,
                          flex: 1.0,
                          sortable: true,
                        }, {
                          id: 'pos_y',
                          header: 'Pos Y',
                          dataIndex: 'pos_y',
                          editor: false,
                          flex: 1.0,
                          sortable: true,
                        }, {
                          id: 'pos_z',
                          header: 'Pos Z',
                          dataIndex: 'pos_z',
                          editor: false,
                          flex: 1.0,
                          sortable: true,
                        },
                      ],
                      plugins: [
                        {ptype: 'cellediting',
                         clicksToEdit: 1},
                      ],
                    }, {
                      xtype: 'button',
                      text: 'Add Keyframe',
                      itemId: 'addKeyframe',
                    }
                  ]
                }, {
                  xtype: 'panel',
                  title: 'Data Editor',
                  itemId: 'data_edit',
                  layout: 'absolute',
                  flex: 1,
                  items : [
                  ]
                },
                ] } /* tabpanel items and entry */
                ]
        }
    ],
});

