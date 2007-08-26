"""
Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

class ReasonInterpreterPanel(wx.Panel):
    def __init__(self, parent, id, locals):
        wx.Panel.__init__(self, parent, id)
        self.shell = wx.py.shell.Shell(
                           parent=self, id=-1, introText="Welcome to Reason.",
                           locals=locals)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.shell, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

def ChooseField(page):
    allFields = page.QueryFields()
    toChoose = nativeFields + [''] + derivedFields
    dlg = wx.SingleChoiceDialog(None,
             'Which field?',
             'Field Chooser (%s)' % outputfile.basename,
             toChoose)
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

