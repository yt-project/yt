/**********************************************************************
The Parameter File Display Widget

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

