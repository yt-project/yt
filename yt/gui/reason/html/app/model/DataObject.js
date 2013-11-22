/**********************************************************************
Data object model

Copyright (c) 2013, yt Development Team.

Distributed under the terms of the Modified BSD License.

The full license is in the file COPYING.txt, distributed with this software.
***********************************************************************/

Ext.define('Reason.model.DataObject', {
    extend: 'Ext.data.Model',
    fields: ['name', 'type', 'filename', 'field_list', 'varname'],
});
