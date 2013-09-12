/**********************************************************************
Main Menu in Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.MainMenu', {
    extend: 'Ext.toolbar.Toolbar',
    alias: 'widget.mainmenu',
    items: [
      {
        text: 'Reason Menu',
        plain: true,
        menu: [
               {xtype: 'menuitem', text: 'Open File', id: 'openfile'},
               {xtype: 'menuitem', text: 'Open Directory', disabled: true},
               {xtype: 'menuseparator'},
               {xtype: 'menuitem', text: 'Save Script',
                id: 'savescript'},
               {xtype: 'menuitem', text: 'Download Script',
                id: 'downloadscript'},
               {xtype: 'menuitem', text: 'Pastebin Script',
                id: 'pastebinscript'},
               {xtype: 'menuseparator'},
               {xtype: 'menuitem', text: 'Enable Debug', id: 'enabledebug'},
               {xtype: 'menuseparator'},
               {xtype:'menuitem', text: 'yt Chat', id: 'ytchat'},
               {xtype: 'menuseparator'},
               {xtype:'menuitem', text: 'Quit', id: 'quit'},
            ],
      },
    ],
});

