/**********************************************************************
File Open Dialog for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.view.MulticastSetup', {
    extend: 'Ext.window.Window',
    alias: 'widget.multicastsetup',
    title: "Set up Multicasting",
    modal: true,
    height: 500,
    width: 800,
    layout: {
        type:'vbox',
        pack: 'start',
        align: 'stretch',
    },
    items: [
        { xtype : "component",
          autoEl : {
              tag : "iframe",
              src : "http://localhost:8080/CreateSession"
          },
          flex: 1.0,
          width: "100%",
        }, {
          xtype: "form",
          labelWidth: 80,
          frame: true,
          height: 120,
          width: "100%",
          items: [
            {
              xtype: 'textfield',
              fieldLabel: 'Session ID',
              itemId: 'session_id',
              width: 300,
            }, {
              xtype: 'textfield',
              fieldLabel: 'Session Token',
              itemId: 'session_token',
              width: 300,
            }, 
          ],
          buttons: [
              {
                  text: 'Multicast', itemId: 'multicast',
              },{
                  text: 'Cancel', itemId: 'cancel',
              }
          ],
        }
    ],
});


