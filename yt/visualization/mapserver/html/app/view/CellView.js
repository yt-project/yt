/**********************************************************************
Notebook Cell view for Reason

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

var cellDisplay = new Ext.XTemplate(
    '<b>Input:</b><br/><br/>',
    '{input}',
    '<br/><br/>',
    '<hr>',
    '<b>Output:</b><br/><br/>',
    '<pre>{output}</pre><br/><br/>'
);

var imageDisplay = new Ext.XTemplate(
    '<b>Image</b><br/><br/>',
    '<hr>',
    '<img src="data:image/png;base64,{image_data}">',
    '<br/><br/>'
);

Ext.define('Reason.view.CellView', {
    extend: 'Ext.grid.Panel',
    alias: 'widget.notebookcells',
    title: 'Cells',
    store: 'CellValues',
    itemId: 'cells',
    autoscroll: true,
    flex: 0.7,
    columns: [{header:'Execution Time', dataIndex: 'executiontime', flex:1}],
    viewConfig: {
        stripeRows: false,
        disableSelection: true,
        trackOver: false,
        enableTextSelection: true,
    },
    features: [{
        ftype: 'rowbody',
        getAdditionalData: function(data, rowIndex, record, orig) {
            var disp;
            console.log(data);
            if(data['image_data'] != '') {
                disp = imageDisplay.apply(data);
            } else {
                disp = cellDisplay.apply(data);
            }
            return {
                rowBody: disp,
                rowBodyCls: 'codeview',
                rowBodyColspan: this.view.headerCt.getColumnCount(),
            };
        }
      }, {
        ftype: 'rowwrap'
      }],
});

