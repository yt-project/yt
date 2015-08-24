/**********************************************************************
Sample widget class

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.controller.widgets.SampleWidget', {
    extend: 'Reason.controller.widgets.BaseWidget',

    statics: {
        widgetName: 'sample_widget',
        supportsDataObjects: false,
        supportsParameterFiles: false,
        displayName: 'Sample Widget',
        preCreation: function(obj) {
            return new this({varname:obj.varname});
        },
    },

    applyPayload: function(payload) {
        return;
    }

});
