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


import yt.lagos as lagos
import yt.raven as raven
import yt.enki as enki
import yt.fido as fido
import yt

import os

#from yt.reason import *
import wx

class OperationsButtonBar(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.parent = parent
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Set up the buttons
        self.sliceButton = wx.Button(self, label = "Slice")
        self.Bind(wx.EVT_BUTTON, self.parent.AddSlice, self.sliceButton)
        self.sizer.Add(self.sliceButton, 1, wx.EXPAND)

        self.projButton = wx.Button(self, label = "Project")
        self.Bind(wx.EVT_BUTTON, self.parent.AddProj, self.projButton)
        self.sizer.Add(self.projButton, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

class ReasonLimitInput(wx.Dialog):
    def __init__(self, plot):
        wx.Dialog.__init__(self, None, -1, 'Limit Setter',
                           size=wx.Size(300,300))
        currentMin, currentMax = plot.norm.vmin, plot.norm.vmax

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.maxSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.maxSizer.AddSpacer(10)
        self.maxSizer.Add(wx.StaticText(self, -1, "Max:"))
        self.maxSizer.AddSpacer(10)
        self.maxVal = wx.TextCtrl(self, -1, "%0.3e" % (currentMax))
        self.maxSizer.Add(self.maxVal)

        self.minSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.minSizer.AddSpacer(10)
        self.minSizer.Add(wx.StaticText(self, -1, "Min:"))
        self.minSizer.AddSpacer(10)
        self.minVal = wx.TextCtrl(self, -1, "%0.3e" % (currentMin))
        self.minSizer.Add(self.minVal)

        self.ok = wx.Button(self, wx.ID_OK, "OK")
        self.ok.SetDefault()

        self.sizer.AddSpacer(10)
        self.sizer.Add(self.maxSizer, 0)
        self.sizer.AddSpacer(10)
        self.sizer.Add(self.minSizer, 0)
        self.sizer.AddSpacer(10)
        self.okSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.okSizer.Add(self.ok, 1, wx.ALIGN_CENTER)
        self.sizer.Add(self.okSizer, 1, wx.ALIGN_CENTER)
        self.sizer.AddSpacer(10)


        self.SetSizer(self.sizer)
        self.Fit()

    def GetData(self):
        return float(self.minVal.GetValue()), \
               float(self.maxVal.GetValue())

class ReasonWidthSelectionWindow(wx.Dialog):
    def __init__(self, outputfile):
        wx.Dialog.__init__(self, None, -1, 'Width Selector',
                           size=wx.Size(300,300))
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.width = wx.TextCtrl(self, -1, "1")
        self.width.SetInsertionPoint(0)
        self.sizer.Add(self.width, 1)

        self.choices = outputfile.units.keys()
        self.choices.sort()
        self.units = wx.Choice(self, -1, (85,25), choices=self.choices)
        self.sizer.Add(self.units, 1)

        self.ok = wx.Button(self, wx.ID_OK, "OK")
        self.sizer.Add(self.ok, 1)
        self.ok.SetDefault()

        self.SetSizer(self.sizer)
        self.Fit()

    def GetData(self):
        return float(self.width.GetValue()), \
               self.choices[self.units.GetSelection()]

def ChooseWidth(outputfile):
    dlg = ReasonWidthSelectionWindow(outputfile)
    resp = dlg.ShowModal()
    w, u = dlg.GetData()
    dlg.Destroy()
    return w, u


