/**********************************************************************
The Plot Window Widget

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
