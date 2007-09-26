"""
Standalone windows and other bits that don't fit elsewhere.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
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

class ReasonParameterFileViewer(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        kwds["title"] = "yt - Reason"
        kwds["size"] = (800,800)
        pf = kwds.pop("outputfile")
        wx.Frame.__init__(self, *args, **kwds)

        # Add the text ctrl
        self.pf = wx.TextCtrl(self, -1, style=wx.TE_READONLY | wx.TE_MULTILINE | wx.HSCROLL)
        self.pf.LoadFile(pf.parameterFilename)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.pf, 1, wx.EXPAND, 0)
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)
        self.Layout()
        self.SetSize((600,600))

class FieldFunctionInspector(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        kwds["title"] = "Field Function Browser"
        kwds["size"] = (600,600)
        wx.Frame.__init__(self, *args, **kwds)

        self.SetupControllers()
        self.DoLayout()
        self.PopulateNotebook()

    def SetupControllers(self):
        self.FunctionNotebook = wx.Choicebook(self, -1)

    def DoLayout(self):
        self.MainSizer = wx.BoxSizer(wx.VERTICAL)

        self.SubMainSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SubMainSizer.AddSpacer(20)
        self.SubMainSizer.Add(self.FunctionNotebook, 1, wx.EXPAND)
        self.SubMainSizer.AddSpacer(20)

        self.MainSizer.AddSpacer(20)
        self.MainSizer.Add(self.SubMainSizer)
        self.MainSizer.AddSpacer(20)

        self.SetSizer(self.MainSizer)
        self.Layout()

    def CreateFieldPage(self, field):
        page = FunctionInspectorPage(self.FunctionNotebook, -1, field=field)
        return page

    def PopulateNotebook(self):
        fs = lagos.fieldInfo.keys()
        fs.sort()
        for field in fs:
            # Check if we want to make a page...
            if lagos.fieldInfo[field][3] == None: continue
            page = self.CreateFieldPage(field)
            self.FunctionNotebook.AddPage(page, field)

class FunctionInspectorPage(wx.Panel):
    def __init__(self, *args, **kwargs):
        if kwargs.has_key("field"): self.field = kwargs.pop("field")
        wx.Panel.__init__(self, *args, **kwargs)

        self.PopulateFields()
        self.SetupControllers()
        self.DoLayout()

    def PopulateFields(self):
        self.units, self.punits, self.logged, self.func = \
            lagos.fieldInfo[self.field]
        self.func_source = lagos.getCode(self.field)

    def SetupControllers(self):
        self.legendUnits = wx.StaticText(self, -1, "Units:", style=wx.ALIGN_RIGHT)
        self.legendProjectedUnits = wx.StaticText(self, -1, "Projected Units:", style=wx.ALIGN_RIGHT)
        self.legendLogged = wx.StaticText(self, -1, "Take Log:", style=wx.ALIGN_RIGHT)

        self.textUnits = wx.StaticText(self, -1, str(self.units), style=wx.ALIGN_LEFT)
        self.textProjectedUnits = wx.StaticText(self, -1, str(self.punits), style=wx.ALIGN_LEFT)

        if self.logged: tl = "Yes"
        else: tl = "No"
        self.textLogged = wx.StaticText(self, -1, tl, style=wx.ALIGN_LEFT)

        self.textFunction = wx.StaticText(self, -1, self.func_source)
        f = self.textFunction.GetFont()
        f.SetFamily(wx.MODERN)
        f = wx.Font(10, wx.MODERN, wx.NORMAL, wx.NORMAL)
        self.textFunction.SetFont(f)

    def DoLayout(self):
        self.MainSizer = wx.BoxSizer(wx.VERTICAL)
        self.UnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ProjectedUnitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.LogSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.FuncSizer = wx.BoxSizer(wx.HORIZONTAL)

        self.UnitSizer.Add(self.legendUnits, 1, wx.EXPAND)
        self.UnitSizer.AddSpacer(30)
        self.UnitSizer.Add(self.textUnits, 1, wx.EXPAND)

        self.ProjectedUnitSizer.Add(self.legendProjectedUnits, 1, wx.EXPAND)
        self.ProjectedUnitSizer.AddSpacer(30)
        self.ProjectedUnitSizer.Add(self.textProjectedUnits, 1, wx.EXPAND)

        self.LogSizer.Add(self.legendLogged, 1, wx.EXPAND)
        self.LogSizer.AddSpacer(30)
        self.LogSizer.Add(self.textLogged, 1, wx.EXPAND)

        self.FuncSizer.Add(self.textFunction, 1, wx.EXPAND)

        self.MainSizer.AddSpacer(10)
        self.MainSizer.Add(self.UnitSizer, 0)
        self.MainSizer.AddSpacer(10)
        self.MainSizer.Add(self.ProjectedUnitSizer, 0)
        self.MainSizer.AddSpacer(10)
        self.MainSizer.Add(self.LogSizer, 0)
        self.MainSizer.AddSpacer(30)
        self.MainSizer.Add(self.FuncSizer, 1, wx.EXPAND)
        self.MainSizer.AddSpacer(10)

        self.SetSizer(self.MainSizer)
        self.Layout()

def ConvertToVTKFloatArray(toConvert):
    nComp = 1
    print "Converting",toConvert.shape
    if len(toConvert.shape) > 1: nComp = toConvert.shape[1]
    nVals = toConvert.shape[0]
    stringValue = toConvert.astype('float32').tostring()
    floatArray = vtk.vtkFloatArray()
    floatArray.SetNumberOfComponents(nComp)
    floatArray.SetVoidArray(stringValue, nVals*nComp, 1)
    floatArrayToReturn = vtk.vtkFloatArray()
    floatArrayToReturn.DeepCopy(floatArray)
    return floatArrayToReturn