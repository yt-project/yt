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

class Profile2DAddField(wx.Dialog):
    def __init__(self, data_object, parent):
        wx.Dialog.__init__(self, parent, -1, title="Add Field to 2D Profile")

        fields = [field for field in QueryFields(data_object)
                  if field not in [parent.data.x_bin_field,
                        parent.data.y_bin_field, parent.data.keys()]]
        
        border = wx.BoxSizer(wx.VERTICAL)
        inner_border = wx.BoxSizer(wx.VERTICAL)

        sbox = wx.StaticBox(self, -1, "Z Field Specifications")
        box = wx.StaticBoxSizer(sbox, wx.VERTICAL)
        gbs = wx.GridBagSizer(5, 5)
        self.z_field = wx.Choice(self, -1, choices=fields, name="Y")
        text = wx.StaticText(self, -1, "Weighting Field")
        self.z_weight = wx.Choice(self, -1, choices=[''] + fields, name="Weight")
        if "CellMassMsun" in fields: self.z_weight.Select(fields.index("CellMassMsun")+1)
        self.z_accx = wx.CheckBox(self, -1, "X Accumulation")
        self.z_accx.SetValue(False)
        self.z_accy = wx.CheckBox(self, -1, "Y Accumulation")
        self.z_accy.SetValue(False)

        gbs.Add(self.z_field, (0,0), (1,2))
        gbs.Add(text, (1,0))
        gbs.Add(self.z_weight, (2,0), (1,2))
        gbs.Add(self.z_accx, (3,0), (1,2))
        gbs.Add(self.z_accy, (4,0), (1,2))
        box.Add(gbs, 1, wx.EXPAND | wx.ALL)
        inner_border.Add(box, 1, wx.EXPAND)
        inner_border.AddSpacer(15)

        gbs = wx.GridBagSizer(5,5)
        ok_button = wx.Button(self, wx.ID_OK, "OK")
        ok_button.SetDefault()
        cancel_button = wx.Button(self, wx.ID_CANCEL, "Cancel")
        gbs.Add(ok_button, (1,0), flag=wx.EXPAND)
        gbs.Add(cancel_button, (1,1), flag=wx.EXPAND)
        inner_border.Add(gbs, 0, wx.EXPAND)
        
        border.Add(inner_border, 1, wx.EXPAND|wx.ALL, 25)
        self.SetSizer(border)
        self.Fit()

    def return_argdict(self):
        argdict = {}
        try:
            argdict['fields'] = str(self.z_field.GetStringSelection())
            argdict['weight'] = str(self.z_weight.GetStringSelection())
            argdict['accumulation'] = (self.z_accx.GetValue(),
                                       self.z_accy.GetValue())
            if argdict['weight'] == '': argdict['weight'] = None
        except ValueError:  
            return None
        return argdict

class Profile2DSetup(wx.Dialog):
    def __init__(self, data_object, parent):
        wx.Dialog.__init__(self, parent, -1, title="Setup 2D Profile")

        fields = QueryFields(data_object)
        
        border = wx.BoxSizer(wx.VERTICAL)
        inner_border = wx.BoxSizer(wx.VERTICAL)

        sbox = wx.StaticBox(self, -1, "X Field Specifications")
        box = wx.StaticBoxSizer(sbox, wx.VERTICAL)
        gbs = wx.GridBagSizer(5, 5)
        self.x_field = wx.Choice(self, -1, choices=fields, name="X")
        text = wx.StaticText(self, -1, "Bin Count")
        self.x_bins = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.x_bins.SetValue('64')
        self.x_bc = wx.CheckBox(self, -1, "Auto-bounds")
        self.x_bc.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.OnToggleXBounds, self.x_bc)
        self.x_min = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.x_max = wx.TextCtrl(self, -1, style=wx.TE_PROCESS_ENTER)
        self.x_log = wx.CheckBox(self, -1, "Take Log")
        self.x_log.SetValue(True)

        gbs.Add(self.x_field, (0,0), (1,2))
        gbs.Add(text, (1,0))
        gbs.Add(self.x_bins, (1,1))
        gbs.Add(self.x_bc, (2,0))
        gbs.Add(self.x_min, (2,1))
        gbs.Add(self.x_max, (3,1))
        gbs.Add(self.x_log, (4,0))
        box.Add(gbs, 1, wx.EXPAND | wx.ALL)
        inner_border.Add(box, 1, wx.EXPAND)
        inner_border.AddSpacer(15)
        
        sbox = wx.StaticBox(self, -1, "Y Field Specifications")
        box = wx.StaticBoxSizer(sbox, wx.VERTICAL)
        gbs = wx.GridBagSizer(5, 5)
        self.y_field = wx.Choice(self, -1, choices=fields, name="Y")
        text = wx.StaticText(self, -1, "Bin Count")
        self.y_bins = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.y_bins.SetValue('64')
        self.y_bc = wx.CheckBox(self, -1, "Auto-bounds")
        self.y_bc.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.OnToggleYBounds, self.y_bc)
        self.y_min = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.y_max = wx.TextCtrl(self, -1, style=wx.TE_PROCESS_ENTER)
        self.y_log = wx.CheckBox(self, -1, "Take Log")
        self.y_log.SetValue(True)

        gbs.Add(self.y_field, (0,0), (1,2))
        gbs.Add(text, (1,0))
        gbs.Add(self.y_bins, (1,1))
        gbs.Add(self.y_bc, (2,0))
        gbs.Add(self.y_min, (2,1))
        gbs.Add(self.y_max, (3,1))
        gbs.Add(self.y_log, (4,0))
        box.Add(gbs, 1, wx.EXPAND | wx.ALL)
        inner_border.Add(box, 1, wx.EXPAND)
        inner_border.AddSpacer(15)

        gbs = wx.GridBagSizer(5,5)
        self.lazy_reader = wx.CheckBox(self, -1, "Memory Conservative")
        self.Bind(wx.EVT_CHECKBOX, self.OnToggleLaziness, self.lazy_reader)
        ok_button = wx.Button(self, wx.ID_OK, "OK")
        ok_button.SetDefault()
        cancel_button = wx.Button(self, wx.ID_CANCEL, "Cancel")
        gbs.Add(self.lazy_reader, (0,0), (1,2), flag=wx.ALIGN_LEFT)
        gbs.Add(ok_button, (1,0), flag=wx.EXPAND)
        gbs.Add(cancel_button, (1,1), flag=wx.EXPAND)
        inner_border.Add(gbs, 0, wx.EXPAND)
        
        self.OnToggleXBounds(None)
        self.OnToggleYBounds(None)

        border.Add(inner_border, 1, wx.EXPAND|wx.ALL, 25)
        self.SetSizer(border)
        self.Fit()

    def __toggle(self, b,mi,ma):
        en = b.GetValue()
        mi.Enable(not en)
        ma.Enable(not en)
        if en:
            mi.SetValue("Min")
            ma.SetValue("Max")
        else:
            mi.SetValue("")
            ma.SetValue("")

    def OnToggleLaziness(self, event):
        if self.lazy_reader.GetValue():
            if self.x_bc.GetValue():
                self.x_bc.SetValue(False)
                self.OnToggleXBounds(None)
            if self.y_bc.GetValue():
                self.y_bc.SetValue(False)
                self.OnToggleYBounds(None)

    def OnToggleXBounds(self, event):
        self.__toggle(self.x_bc, self.x_min, self.x_max)

    def OnToggleYBounds(self, event):
        self.__toggle(self.y_bc, self.y_min, self.y_max)

    def return_argdict(self):
        argdict = {}
        try:
            if self.lazy_reader.GetValue():
                argdict['lazy_reader'] = True
                argdict['x_lower_bound'] = float(self.x_min.GetValue())
                argdict['x_upper_bound'] = float(self.x_max.GetValue())
                argdict['y_lower_bound'] = float(self.y_min.GetValue())
                argdict['y_upper_bound'] = float(self.y_max.GetValue())
            else:
                argdict['lazy_reader'] = False
                argdict['x_lower_bound'] = None
                argdict['x_upper_bound'] = None
                argdict['y_lower_bound'] = None
                argdict['y_upper_bound'] = None
            argdict['x_log'] = self.x_log.GetValue()
            argdict['y_log'] = self.y_log.GetValue()
            argdict['x_n_bins'] = int(self.x_bins.GetValue())
            argdict['y_n_bins'] = int(self.y_bins.GetValue())
            argdict['x_bin_field'] = str(self.x_field.GetStringSelection())
            argdict['y_bin_field'] = str(self.y_field.GetStringSelection())
        except ValueError:  
            return None
        return argdict

class ProjectionSetup(wx.Dialog):
    def __init__(self, data_object, parent):
        wx.Dialog.__init__(self, parent, -1, title="Setup Projection")

        fields = QueryFields(data_object, True)
        
        border = wx.BoxSizer(wx.VERTICAL)
        gbs = wx.GridBagSizer(5, 5)
        text = wx.StaticText(self, -1, "Projected Field")
        self.field = wx.Choice(self, -1, choices=fields, name="Field")
        if "Density" in fields: self.field.Select(fields.index("Density"))
        gbs.Add(text, (0,0))
        gbs.Add(self.field, (0,1))

        text = wx.StaticText(self, -1, "Weighting Field")
        self.weight_field = wx.Choice(self, -1, choices=[""]+fields,
                                      name="Weight field")
        gbs.Add(text, (1,0))
        gbs.Add(self.weight_field, (1,1))

        gbs.Add(wx.StaticText(self, -1, "Axes to Project"), (2,0),
                flag = wx.ALIGN_CENTER | wx.ALL)
        box = wx.BoxSizer(wx.HORIZONTAL)
        self.x_ax = wx.CheckBox(self, -1, "X")
        self.y_ax = wx.CheckBox(self, -1, "Y")
        self.z_ax = wx.CheckBox(self, -1, "Z")
        for i in [self.x_ax, self.y_ax, self.z_ax]:
            box.Add(i, 1, wx.EXPAND, border=5)
            i.SetValue(True)
        gbs.Add(box, (2,1), flag = wx.ALL | wx.EXPAND)

        box = wx.BoxSizer(wx.HORIZONTAL)
        ok_button = wx.Button(self, wx.ID_OK, "OK")
        ok_button.SetDefault()
        cancel_button = wx.Button(self, wx.ID_CANCEL, "Cancel")
        box.Add(ok_button, 1, wx.EXPAND)
        box.Add(cancel_button, 1, wx.EXPAND)
        gbs.Add(box, (3,1), flag=wx.EXPAND)

        border.Add(gbs, 1, wx.EXPAND|wx.ALL, 25)
        self.SetSizer(border)
        self.Fit()
        


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
        field = lagos.fieldInfo[self.field]
        self.units = field.get_units()
        self.punits = field.get_projected_units()
        self.logged = field.take_log
        self.func_source = lagos.fieldInfo[self.field].get_source()

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


class ReasonEditorNotebookFrame(wx.py.editor.EditorNotebookFrame):

    def _setup(self):
        """Setup prior to first buffer creation.
        Called automatically by base class during init.

        Mostly taken from the wx.py.editor source.
        """
        self.notebook = wx.py.editor.EditorNotebook(parent=self)
        intro = 'Reason 0.3'
        if hasattr(self.Parent, 'locals'):
            namespace = self.Parent.locals
        else:
            import imp
            module = imp.new_module('__main__')
            import __builtin__
            module.__dict__['__builtins__'] = __builtin__
            namespace = module.__dict__.copy()
        self.crust = wx.py.crust.Crust(parent=self.notebook, intro=intro, locals=namespace)
        self.shell = self.crust.shell
        # Override the filling so that status messages go to the status bar.
        self.crust.filling.tree.setStatusText = self.SetStatusText
        # Override the shell so that status messages go to the status bar.
        self.shell.setStatusText = self.SetStatusText
        # Fix a problem with the sash shrinking to nothing.
        self.crust.filling.SetSashPosition(40)
        self.notebook.AddPage(page=self.crust, text='*Shell*', select=True)
        self.setEditor(self.crust.editor)
        self.crust.editor.SetFocus()
        self.MenuBar.GetMenu(0).GetMenuItems()[11].SetText("Run Script")
