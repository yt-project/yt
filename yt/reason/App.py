"""
Main application for Reason.  Includes the basic window outline.

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

_StaticOutputMenuItems = ["proj","slice"]
_SphereObjectMenuItems = ["phase"]
_ProjObjectMenuItems = []
_SliceObjectMenuItems = []

class ReasonMainWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        kwds["title"] = "yt - Reason"
        kwds["size"] = (ytcfg.getint("reason","width"),
                        ytcfg.getint("reason","height"))
        wx.Frame.__init__(self, *args, **kwds)

        self.__setup_controls()

        self.__setup_menubar()
        self.__setup_popup_menu()
        self.__setup_toolbar()
        self.__setup_data_tree()

        self.status_bar = self.CreateStatusBar(4, 0)

        self.__set_properties()
        self.__do_layout()

        Publisher().subscribe(self.MessagePageDeleted, ('page_deleted'))

    def __setup_controls(self):

        self.main_splitter = wx.SplitterWindow(self, -1)
        self.view_panel = wx.SplitterWindow(self.main_splitter, -1)
        self.data_panel = wx.Panel(self.view_panel, -1)
        self.data_tree = wx.TreeCtrl(self.data_panel, -1,
                style=wx.TR_HIDE_ROOT | wx.TR_LINES_AT_ROOT | wx.TR_HAS_BUTTONS)
        self.plot_panel = PlotPanel(parent=self.view_panel)
        self.__setup_interpreter()

        self.main_splitter.SetMinimumPaneSize(20)
        self.main_splitter.SplitHorizontally(self.view_panel, self.int_panel, -100)
        self.view_panel.SetMinimumPaneSize(20)
        self.view_panel.SplitVertically(self.data_panel, self.plot_panel, 200)

    def __setup_interpreter(self):
        self.windows = []
        self.outputs = []
        self.data_objects = []
        self.locals = {'lagos':lagos,
                       'raven':raven,
                       'enki':enki,
                       'raven':raven,
                       'outputs':self.outputs,
                       'windows':self.windows,
                       'mainwindow':self,
                       'data_objects':self.data_objects}
        self.int_panel = wx.Panel(self.main_splitter, -1)
        self.interpreter = ReasonInterpreterPanel(self.int_panel, -1, self.locals)

    def __setup_menubar(self):
        menu_bar = wx.MenuBar()
        file_menu = wx.Menu()
        menu_bar.Append(file_menu, "File")

        # Set up IDs for event binding

        open_hierarchy = file_menu.Append(-1, "Open Hierarchy")
        field_inspector = file_menu.Append(-1, "Inspect Fields")
        save_image = file_menu.Append(-1, "Save Image")
        file_menu.AppendSeparator()
        exit = file_menu.Append(-1, "Exit")

        self.Bind(wx.EVT_MENU, self.OnOpenHierarchy, open_hierarchy)
        self.Bind(wx.EVT_MENU, self.OnInspectFields, field_inspector)
        self.Bind(wx.EVT_MENU, self.OnSaveImage, save_image)
        self.Bind(wx.EVT_MENU, self.OnExit, exit)

        self.SetMenuBar(menu_bar)

    def __setup_popup_menu(self):
        self.PopupMenu = wx.Menu()
        self.PopupMenuIds = {}
        self.PopupMenuIds["slice"] = self.PopupMenu.Append(-1, "Slice")
        self.PopupMenuIds["proj"] = self.PopupMenu.Append(-1, "Project")
        self.PopupMenuIds["phase"] = self.PopupMenu.Append(-1, "Phase Plot")
        self.PopupMenuIds["vr"] = self.PopupMenu.Append(-1, "Volume Render")

        self.Bind(wx.EVT_MENU, self._add_slice, self.PopupMenuIds["slice"])
        self.Bind(wx.EVT_MENU, self._add_proj, self.PopupMenuIds["proj"])
        self.Bind(wx.EVT_MENU, self._add_phase, self.PopupMenuIds["phase"])

    def __setup_data_tree(self):

        self.root = self.data_tree.AddRoot("You shouldn't see me!")
        self.fido_root = self.data_tree.AppendItem(self.root, "Stored Outputs")
        self.output_root = self.data_tree.AppendItem(self.root, "Selected Outputs")
        self.data_root = self.data_tree.AppendItem(self.root, "Data Objects")

        self.data_tree.Expand(self.fido_root)
        self.data_tree.Expand(self.output_root)
        self.data_tree.Expand(self.data_root)

        self.data_tree.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
        self.Bind(wx.EVT_TREE_ITEM_EXPANDED, self.OnItemExpanded, self.data_tree)

        self.__setup_fido_tree()

    def __setup_fido_tree(self, event=None):
        # Calling this delete may not be wise.
        # However, as long as we have a distinction between
        # the data outputs and the created data objects, it should be okay.
        self.data_tree.DeleteChildren(self.fido_root)
        gc = fido.GrabCollections()
        for c in gc:
            cRoot = self.data_tree.AppendItem(self.fido_root, c.title)
            for fn in c:
                if not os.path.isfile(fn): continue
                try:
                    z = str(fido.get_parameter_line(fn,
                                 "CosmologyCurrentRedshift"))
                except:
                    z = "N/A"
                tt = str(fido.get_parameter_line(fn,
                             "InitialTime"))
                tid = wx.TreeItemData((fn, tt, z, _StaticOutputMenuItems))
                ni = self.data_tree.AppendItem(cRoot,
                    "%s" % (os.path.basename(fn)), data=tid)

    def __setup_toolbar(self):
        # Tool Bar
        self._VMTB_REREADFIDO = wx.NewId()
        self._VMTB_FULLDOMAIN = wx.NewId()
        self._VMTB_CHANGEZOOM = wx.NewId()
        self._VMTB_REDRAW = wx.NewId()
        self._VMTB_SAVE = wx.NewId()
        self._VMTB_FIELDSWITCH = wx.NewId()
        self._VMTB_CHANGELIMITS = wx.NewId()
        self._VMTB_VIEWPF = wx.NewId()
        self._VMTB_VELPLOT = wx.NewId()

        self.toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL|wx.TB_TEXT|wx.TB_HORZ_LAYOUT)
        font = self.toolbar.GetFont()
        font.SetFamily(wx.MODERN)
        self.toolbar.SetFont(font)
        self.toolbar.SetToolBitmapSize((16,16))

        self.SetToolBar(self.toolbar)
        def AddButton(id, label, tooltip="", bitmapID=None):
            if bitmapID != None:
                bm = wx.ArtProvider.GetBitmap(bitmapID, wx.ART_TOOLBAR, (16,16))
            else: bm = wx.NullBitmap
            self.toolbar.AddLabelTool(id, label, bm, bm, wx.ITEM_NORMAL, tooltip, "")
            self.toolbar.AddSeparator()

        self.toolbar.AddSeparator()
        AddButton(self._VMTB_REREADFIDO,"Update OutputList",
                                        "Reread from the Fido database", wx.ART_TIP)
        AddButton(self._VMTB_REDRAW,"Redraw", "Force a redraw", wx.ART_REDO)
        self.available_fields = wx.Choice(self.toolbar, id=self._VMTB_FIELDSWITCH, choices = [])
        self.toolbar.AddControl(self.available_fields)
        Publisher().subscribe(self.MessageUpdateToolbarFields, ('page_changed'))
        AddButton(self._VMTB_FULLDOMAIN, "Zoom Top",  "Zoom to the top level", wx.ART_FIND)
        AddButton(self._VMTB_CHANGELIMITS, "Change Limits", "Change the colorbar limits", wx.ART_GO_UP)
        AddButton(self._VMTB_VIEWPF, "View ParameterFile", "View the parameter file", wx.ART_NORMAL_FILE)
        cl = wx.ArtProvider.GetBitmap(wx.ART_TICK_MARK, wx.ART_TOOLBAR, (16,16))
        #self.toolbar.AddCheckLabelTool(self._VMTB_VELPLOT, "VelVecs", cl, shortHelp="Plot Velocity Vectors")
        self.toolbar.AddSeparator()

        self.Bind(wx.EVT_MENU, self.__setup_fido_tree, id=self._VMTB_REREADFIDO)
        self.Bind(wx.EVT_CHOICE, self.plot_panel.OnCallSwitchField, id=self._VMTB_FIELDSWITCH)
        self.Bind(wx.EVT_MENU, self.plot_panel.OnCallSetWidth, id=self._VMTB_CHANGEZOOM)
        self.Bind(wx.EVT_MENU, self.plot_panel.OnCallRedraw, id=self._VMTB_REDRAW)
        self.Bind(wx.EVT_MENU, self.plot_panel.OnCallZoomTop, id=self._VMTB_FULLDOMAIN)
        self.Bind(wx.EVT_MENU, self.plot_panel.OnCallSetZLim, id=self._VMTB_CHANGELIMITS)
        self.Bind(wx.EVT_MENU, self.plot_panel.OnCallViewPF, id=self._VMTB_VIEWPF)

    def __set_properties(self):
        self.toolbar.SetToolBitmapSize((24, 24))
        self.toolbar.Realize()
        self.status_bar.SetStatusWidths([-1,-1,-1,-1])

    def __do_layout(self):
        MainWindowSizer = wx.BoxSizer(wx.VERTICAL)
        DataPanelSizer = wx.BoxSizer(wx.VERTICAL)

        IntPanelSizer = wx.BoxSizer(wx.HORIZONTAL)
        IntPanelSizer.Add(self.interpreter, 1, wx.EXPAND, 0)
        self.int_panel.SetSizer(IntPanelSizer)

        DataPanelSizer.Add(self.data_tree, 1, wx.EXPAND, 0)
        self.data_panel.SetSizer(DataPanelSizer)
        self.data_panel.Layout()

        MainWindowSizer.Add(self.main_splitter, 1, wx.EXPAND)
        self.SetSizer(MainWindowSizer)

        self.Layout()

    def _update_toolbar_fields(self, page):
        newItems = None
        if hasattr(page,'QueryFields'): newItems = page.QueryFields()
        if not newItems:
            self.available_fields.Enable(False)
        else:
            self.available_fields.Enable(True)
            self.available_fields.SetItems(newItems)

    def _add_static_output(self, filename):
        # Alright, we choose the hierarchy in the file selector,
        # so let's strip that extension off
        fn = filename[:-10]
        eso = lagos.EnzoStaticOutput(fn)
        try:
            z = str(eso["CosmologyCurrentRedshift"])
        except:
            z = "N/A"
        tid = wx.TreeItemData((eso, str(eso["InitialTime"]), z, _StaticOutputMenuItems))
        self.outputs.append(eso)
        ni = self.data_tree.AppendItem(self.output_root, "%s" % (eso.basename), data=tid)
        self.data_tree.Expand(self.output_root)

    def _add_data_object(self, title, object, mids):
        self.data_objects.append(object)
        tid = wx.TreeItemData((object, title, len(self.data_objects), mids))
        ni = self.data_tree.AppendItem(self.data_root, "%s" % (title), data=tid)
        self.data_tree.Expand(self.data_root)

    def _add_sphere(self, title, sphere):
        # These all get passed in
        self._add_data_object(title, sphere, _SphereObjectMenuItems)

    def _add_phase(self, event=None):
        MyID = wx.NewId()
        self.interpreter.shell.writeOut("\n")
        o = self.get_output()
        t = "Phase Plot"
        self.windows.append( \
            PhasePlotPage(parent=self.plot_panel.nb,
                          status_bar=self.status_bar,
                          dataObject = o,
                          CreationID = MyID,
                          mw = self))
        self.interpreter.shell.writeOut("Adding phase plot\n")
        self.plot_panel.AddPlot(self.windows[-1], t, MyID)
        mylog.debug("Adding with ID: %s", MyID)
        self.interpreter.shell.push("\n")

    def _add_proj(self, event=None):
        MyID = wx.NewId()
        self.interpreter.shell.writeOut("\n")
        o = self.get_output()
        field = "Density"
        width = 1.0
        unit = "1"
        for i, ax in zip(range(3), 'xyz'):
            t = "%s - Projection - %s" % (o.basename, ax)
            self.interpreter.shell.writeOut("Adding %s projection of %s\n" % (ax, o))
            self.windows.append( \
                ProjPlotPage(parent=self.plot_panel.nb,
                              status_bar=self.status_bar,
                              outputfile = o,
                              axis=i,
                              field = field,
                              mw = self, CreationID=MyID))
            self.plot_panel.AddPlot(self.windows[-1], t, MyID)
            self._add_data_object("Proj: %s %s" % (o, ax),
                               self.windows[-1].plot.data,
                               _ProjObjectMenuItems)
            print "Adding with ID:", MyID
        for w in self.windows[-3:]: w.ChangeWidth(1,'1')
        self.interpreter.shell.push("\n")

    def _add_slice(self, event=None):
        MyID = wx.NewId()
        self.interpreter.shell.writeOut("\n")
        o = self.get_output()
        field = "Density"
        width = 1.0
        unit = "1"
        for i, ax in zip(range(3), 'xyz'):
            t = "%s - Slice - %s" % (o.basename, ax)
            self.interpreter.shell.writeOut("Adding %s slice of %s\n" % (ax, o))
            self.windows.append( \
                SlicePlotPage(parent=self.plot_panel.nb,
                              status_bar=self.status_bar,
                              outputfile = o,
                              axis=i,
                              field = field,
                              mw = self, CreationID=MyID))
            self.plot_panel.AddPlot(self.windows[-1], t, MyID)
            self._add_data_object("Slice: %s %s" % (o, ax),
                               self.windows[-1].plot.data,
                               _SliceObjectMenuItems)
            print "Adding with ID:", MyID
        for w in self.windows[-3:]: w.ChangeWidth(1,'1')
        self.interpreter.shell.push("\n")

    def get_output(self, event=None):
        # Figure out which outputs are selected
        #tid = self.data_tree.GetFirstSelected()
        tid = self.data_tree.GetSelection()
        ii = self.data_tree.GetItemData(tid).GetData()[0]
        if isinstance(ii, types.StringTypes):
            ii = lagos.EnzoStaticOutput(ii) # Instantiate here
            self.outputs.append(ii)
            fn, z, t, mids = self.data_tree.GetItemData(tid).GetData()
            newData = wx.TreeItemData((ii, z, t, mids))
            self.data_tree.SetItemData(tid, newData)
        print "Got output:", ii
        return ii

    # Functions bound to messages in pubsub

    def MessageUpdateToolbarFields(self, message):
        page = message.data
        self._update_toolbar_fields(page)

    def MessagePageDeleted(self, message):
        id = message.data
        del self.windows[id]

    # Functions bound exclusively to events:

    def OnExit(self, event):
        self.Close()

    def OnInspectFields(self, event):
        p = FieldFunctionInspector(self, -1)
        p.Show()

    def OnOpenHierarchy(self, event):
        wildcard = "Hierarchy (*.hierarchy)|*.hierarchy|" \
                   "All files (*,*)|*.*"
        dialog = wx.FileDialog(None, "Choose your hierarchy", os.getcwd(),
                               "", wildcard, wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            file = dialog.GetPath()
            print file
            self._add_static_output(file)
            #self.RefreshOutputs()
        dialog.Destroy()

    def OnSaveImage(self, event):
        pgI = self.plot_panel.nb.Selection
        pg = self.plot_panel.nb.GetPage(pgI)
        if not hasattr(pg, 'SaveImage'): return
        pg.SaveImage()

    def OnItemExpanded(self, event):
        if event.GetItem() == self.fido_root:
            mylog.info("Reloading fido outputs")
            self.__setup_fido_tree()

    def OnRightDown(self, event):
        pt = event.GetPosition();
        item, flags = self.data_tree.HitTest(pt)
        if item:
            self.data_tree.SelectItem(item)
            pos = event.GetPosition()
            self.ContextMenuPosition = pos
            itemData = self.data_tree.GetItemData(item).Data
            if not itemData: return
            for n,d in self.PopupMenuIds.items():
                self.PopupMenu.Enable(d.Id,False)
            if itemData:
                for n in itemData[3]:
                    self.PopupMenu.Enable(self.PopupMenuIds[n].Id, True)
            #self.PopupMenu.Enable(self.PopupMenuIds["proj"].Id,False)
            self.data_tree.PopupMenu(self.PopupMenu, pos)
