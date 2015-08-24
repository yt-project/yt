/**********************************************************************
The Plot Window Widget

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define("Reason.controller.widgets.ProgressBar", {
    extend: 'Reason.controller.widgets.BaseWidget',
    requires: [],

    templates: {

    },

    createView: function() {
        Ext.MessageBox.show({
            title: 'yt is working ...',
            msg: this.payload['data'].title,
            progressText: 'Progress',
            width: 300,
            progress: true,
            closable: false,
        });
    },

    applyPayload: function(payload) {
        var i = payload['value'];
        if (i == -1) {
            Ext.MessageBox.hide();
        } else {
            Ext.MessageBox.updateProgress(i, Math.round(100*i)+'% completed');
        }
    },

    statics: {
        widgetName: "progressbar",
        supportsDataObjects: false,
        supportsParameterFiles: false,
        displayName: 'Do not use',
    },

});
