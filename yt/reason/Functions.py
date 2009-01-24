"""
Standalone functions which may or may not create and interact with windows,
but not any subclassing of window objects.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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
"""

from yt.reason import *

subscription_pairs = [("width","ChangeWidthFromMessage"),
                      ("field","ChangeFieldFromMessage"),
                      ("limits","ChangeLimitsFromMessage"),
                      ("center","ChangeCenterFromMessage"),
                      ("field_param","ChangeFieldParamFromMessage"),
                      ("wipe","WipePlotDataFromMessage"),
                      ("cmap","ChangeColorMapFromMessage")]

def QueryFields(outputfile, only_display_fields = False):
    fields = []
    for f in outputfile.hierarchy.derived_field_list:
        if f in outputfile.field_info and \
            outputfile.field_info[f].particle_type: continue
        if only_display_fields and \
            f in outputfile.field_info and \
            not outputfile.field_info[f].display_field: continue
        fields.append(f)
    return sorted(fields)

def ChooseField(data_object):
    fields = QueryFields(data_object)
    dlg = wx.SingleChoiceDialog(None,
             'Which field?',
             'Field Chooser (%s)' % data_object,
             fields)
    response = None
    if dlg.ShowModal() == wx.ID_OK:
        response = dlg.GetStringSelection()
    if response == "":
        response = None
    return response

def ChooseLimits(plot):
    dlg = ReasonLimitInput(plot)
    resp = dlg.ShowModal()
    zmin, zmax = dlg.GetData()
    dlg.Destroy()
    return zmin, zmax

def get_new_updateNamespace(my_locals):
    def updateNamespace(self):
        """Update the namespace for autocompletion and calltips.

        Return True if updated, False if there was an error."""
        if not self.interp or not hasattr(self.editor, 'getText'):
            return False
        syspath = sys.path
        sys.path = self.syspath
        text = self.editor.getText()
        text = text.replace('\r\n', '\n')
        text = text.replace('\r', '\n')
        name = self.modulename or self.name
        try:
            try:
                code = compile(text, name, 'exec')
            except:
                raise
            try:
                exec code in my_locals
            except:
                raise
        finally:
            sys.path = syspath
            for m in sys.modules.keys():
                if m not in self.modules:
                    del sys.modules[m]
    return updateNamespace

