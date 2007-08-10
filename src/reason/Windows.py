import yt.lagos as lagos
import yt.raven as raven
import yt.enki as enki
import yt.fido as fido
import yt
from yt.arraytypes import *

from math import log10, sqrt

import os, Toolbars

import yt.raven.backends.MPL as be
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx, wx.py
import matplotlib.backends.backend_wx as be_wx

class ReasonMainWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        kwds["title"] = "yt - Reason"
        kwds["size"] = (800,800)
        wx.Frame.__init__(self, *args, **kwds)

        self.windows = []
        self.outputs = []
        self.locals = {'lagos':lagos,
                       'raven':raven,
                       'enki':enki,
                       'raven':raven,
                       'outputs':self.outputs,
                       'windows':self.windows}

        self.outputPanel = wx.Panel(self, -1)
        self.intPanel = wx.Panel(self, -1)
        self.interpreter = ReasonInterpreterPanel(self.intPanel, -1, self.locals)
        self.outputList = wx.ListCtrl(self.outputPanel, -1, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        self.SliceButton = wx.Button(self.outputPanel, -1, "Slice")
        self.ProjectButton = wx.Button(self.outputPanel, -1, "Project")
        self.button_2 = wx.Button(self.outputPanel, -1, "button_2")
        self.button_3 = wx.Button(self.outputPanel, -1, "button_3")

        self.SetupMenubar()

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: ReasonMainWindow.__set_properties
        self.outputList.SetMinSize((300,300))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ReasonMainWindow.__do_layout
        MainWindowSizer = wx.GridSizer(1, 2, 0, 0)
        OutputPanelSizer = wx.BoxSizer(wx.VERTICAL)
        ButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        IntPanelSizer = wx.BoxSizer(wx.HORIZONTAL)
        IntPanelSizer.Add(self.interpreter, 1, wx.EXPAND, 0)
        self.intPanel.SetSizer(IntPanelSizer)
        MainWindowSizer.Add(self.intPanel, 1, wx.EXPAND, 0)
        OutputPanelSizer.Add(self.outputList, 1, wx.EXPAND, 0)
        ButtonSizer.Add(self.SliceButton, 0, 0, 0)
        ButtonSizer.Add(self.ProjectButton, 0, 0, 0)
        ButtonSizer.Add(self.button_2, 0, 0, 0)
        ButtonSizer.Add(self.button_3, 0, 0, 0)
        OutputPanelSizer.Add(ButtonSizer, 0, 0, 0)
        self.outputPanel.SetSizer(OutputPanelSizer)
        MainWindowSizer.Add(self.outputPanel, 1, wx.EXPAND, 0)
        self.SetSizer(MainWindowSizer)
        MainWindowSizer.Fit(self)
        self.Layout()

        self.Bind(wx.EVT_BUTTON, self.AddSlice, self.SliceButton)
        self.Bind(wx.EVT_BUTTON, self.AddProj, self.ProjectButton)

        #self.AddStaticOutputFile("/Users/matthewturk/Research/data/mornkr/galaxy0398.dir/galaxy0398.hierarchy")
        self.AddStaticOutputFile("/c/britton/data/cl-4/DataDir0036/DataDump0036.hierarchy")
        #self.AddStaticOutputFile("/c/mturk/data/RD0069/RD0069.hierarchy")
        self.RefreshOutputs()

    def SetupMenubar(self):
        menuBar = wx.MenuBar()
        fileMenu = wx.Menu()
        menuBar.Append(fileMenu, "File")
        
        # Set up IDs for event binding

        openHierarchy = fileMenu.Append(-1, "Open hierarchy")
        fileMenu.AppendSeparator()
        exit = fileMenu.Append(-1, "Exit")

        self.Bind(wx.EVT_MENU, self.OnOpenHierarchy, openHierarchy)
        self.Bind(wx.EVT_MENU, self.OnExit, exit)

        self.SetMenuBar(menuBar)

    def OnExit(self, event):
        self.Close()

    def OnOpenHierarchy(self, event):
        wildcard = "Hierarchy (*.hierarchy)|*.hierarchy|" \
                   "All files (*,*)|*.*"
        dialog = wx.FileDialog(None, "Choose your hierarchy", os.getcwd(),
                               "", wildcard, wx.OPEN)
        if dialog.ShowModal() == wx.ID_OK:
            file = dialog.GetPath()
            print file
            self.AddStaticOutputFile(file)
            self.RefreshOutputs()
        dialog.Destroy()

    def AddStaticOutputFile(self, filename):
        # Alright, we choose the hierarchy in the file selector,
        # so let's strip that extension off
        fn = filename[:-10]
        eso = lagos.EnzoStaticOutput(fn, data_style=5)
        self.outputs.append(eso)
        # Note that we do NOT refresh here, as we may ultimately want to add
        # several at once.  Refreshing the outputs is left as an exercise to
        # the caller.

    def RefreshOutputs(self, event=None):
        self.outputList.ClearAll()
        self.outputList.InsertColumn(0, "Basename", width=100)
        self.outputList.InsertColumn(1, "Time", width=100)
        self.outputList.InsertColumn(2, "Redshift", width=100)
        self.outputList.InsertColumn(3, "index", width=50)
        for i,o in enumerate(self.outputs):
            self.outputList.InsertStringItem(i, o.basename)
            self.outputList.SetStringItem(i, 1, str(o["InitialTime"]))
            try:
                z = str(o["CosmologyCurrentRedshift"])
            except:
                z = "N/A"
            self.outputList.SetStringItem(i, 2, z)
            self.outputList.SetStringItem(i, 3, str(i))

    def AddProj(self, event=None):
        for o in self.GetOutputs():
            field = Toolbars.ChooseField(o)
            if not field:
                continue
            #width, unit = Toolbars.ChooseWidth(o)
            width = 1.0
            unit = "1"
            for i, ax in zip(range(1), 'xyz'):
                t = "%s - Projection - %s" % (o.basename, ax)
                self.windows.append( \
                    ReasonVMPlotFrame(parent=self, title=t,
                                      axis=i, size=(500,500),
                                      outputfile = o,
                                      vmtype=ReasonProjPlotPanel,
                                      field = field))
                self.windows[-1].Show()
                self.windows[-1].vmPanel.plot.set_width(width, unit)

    def AddSlice(self, event=None):
        for o in self.GetOutputs():
            #field = Toolbars.ChooseField(o)
            field = "Density"
            if not field:
                continue
            #width, unit = Toolbars.ChooseWidth(o)
            width = 1.0
            unit = "1"
            for i, ax in zip(range(3), 'xyz'):
                t = "%s - Slice - %s" % (o.basename, ax)
                self.windows.append( \
                    ReasonVMPlotFrame(parent=self, title=t,
                                      axis=i, size=(500,500),
                                      outputfile = o,
                                      vmtype=ReasonSlicePlotPanel,
                                      field = field))
                self.windows[-1].Show()
                self.windows[-1].vmPanel.plot.set_width(width, unit)

    def GetOutputs(self, event=None):
        # Figure out which outputs are selected
        os = []
        k = self.outputList.GetFirstSelected()
        while k != -1:
            os.append(self.outputs[k])
            k = self.outputList.GetNextSelected(k)
        return os


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

class ReasonVMPlotFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: ReasonVMPlotFrame.__init__
        self.CurrentlyResizing = False
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        if not kwds.has_key("vmtype"):
            vmtype = ReasonSlicePlotPanel
        else:
            vmtype = kwds.pop("vmtype")
        outputfile = kwds.pop("outputfile")
        self.outputfile = outputfile
        field = kwds.pop("field")
        axis = kwds.pop("axis")
        self.parent = kwds["parent"]
        wx.Frame.__init__(self, *args, **kwds)
        self.PlotPanel = wx.Panel(self, -1)

        # Tool Bar
        self._VMTB_FULLDOMAIN = wx.NewId()
        self._VMTB_CHANGEZOOM = wx.NewId()
        self._VMTB_REDRAW = wx.NewId()
        self._VMTB_SAVE = wx.NewId()
        self._VMTB_FIELDSWITCH = wx.NewId()
        self._VMTB_CHANGELIMITS = wx.NewId()
        self._VMTB_VIEWPF = wx.NewId()

        self.toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL|wx.TB_NOICONS|wx.TB_HORZ_LAYOUT)
        font = self.toolbar.GetFont()
        font.SetFamily(wx.MODERN)
        self.toolbar.SetFont(font)

        self.SetToolBar(self.toolbar)
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_REDRAW, "Redraw", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "Force a redraw", "")
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_FIELDSWITCH, "Change Field", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "Change the displayed field", "")
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_CHANGEZOOM, "Change Width", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "Change the displayed width", "")
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_FULLDOMAIN, "Zoom Top", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "Zoom to the top level", "")
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_CHANGELIMITS, "Change Limits", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "Change the colorbar limits", "")
        self.toolbar.AddSeparator()
        self.toolbar.AddLabelTool(self._VMTB_VIEWPF, "View ParameterFile", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "View the parameter file", "")
        self.toolbar.AddSeparator()
        # Tool Bar end
        self.statusBar = self.CreateStatusBar(1, 0)

        self.vmPanel = vmtype(outputfile, self.statusBar, self.PlotPanel, -1, 
                                       field=field, axis=axis)
        

        #self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_MENU, self.vmPanel.switch_field, id=self._VMTB_FIELDSWITCH)
        self.Bind(wx.EVT_MENU, self.vmPanel.set_width, id=self._VMTB_CHANGEZOOM)
        self.Bind(wx.EVT_MENU, self.vmPanel.redraw, id=self._VMTB_REDRAW)
        self.Bind(wx.EVT_MENU, self.vmPanel.fulldomain, id=self._VMTB_FULLDOMAIN)
        self.Bind(wx.EVT_MENU, self.vmPanel.set_zlim, id=self._VMTB_CHANGELIMITS)
        self.Bind(wx.EVT_MENU, self.viewPF, id=self._VMTB_VIEWPF)

        self.__set_properties()
        self.__do_layout()
        self.SetSize((600,600))
        #self.OnSize()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: ReasonVMPlotFrame.__set_properties
        self.toolbar.SetToolBitmapSize((24, 24))
        self.toolbar.Realize()
        self.statusBar.SetStatusWidths([-1])
        # statusbar fields
        statusBar_fields = ["statusBar"]
        for i in range(len(statusBar_fields)):
            self.statusBar.SetStatusText(statusBar_fields[i], i)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ReasonVMPlotFrame.__do_layout
        PlotPanelSizer = wx.BoxSizer(wx.VERTICAL)
        vmSizer = wx.BoxSizer(wx.VERTICAL)
        vmSizer.Add(self.vmPanel, 1, wx.EXPAND, 0)
        self.PlotPanel.SetSizer(vmSizer)
        PlotPanelSizer.Add(self.PlotPanel, 1, wx.EXPAND, 0)
        self.SetSizer(PlotPanelSizer)
        PlotPanelSizer.Fit(self)
        self.Layout()
        # end wxGlade

    def viewPF(self, event=None):
        r = ReasonParameterFileViewer(self, -1, outputfile=self.outputfile)
        r.Show()

    def OnSize(self, event):
        if self.CurrentlyResizing == True:
            return
        self.CurrentlyResizing = True
        #w,h = self.GetEffectiveMinSize()
        w,h = event.GetSize()
        xx = min(w,h)
        self.SetSize((xx,xx))
        w,h = self.vmPanel.GetClientSize()
        xx = min(w,h)
        xy = xx / self.vmPanel.figure.get_dpi()
        self.vmPanel.figure.set_size_inches(xy,xy)
        self.vmPanel.figure_canvas.SetSize((xx,xx))
        self.vmPanel.plot.redraw_image()
        self.Layout()
        self.CurrentlyResizing = False

class ReasonVMPlotPanel(wx.Panel):
    def __init__(self, outputfile, statusBar, parent=None, title="VM Plot", axis=0, size=(800,800),
                 field = "Density"):
        wx.Panel.__init__(self, parent, title, size)

        self.AmDrawingCircle = False

        self.outputfile = outputfile
        self.figure = be.matplotlib.figure.Figure((4,4))
        self.axes = self.figure.add_subplot(111)
        self.statusBar = statusBar

        self.field = field

        self.SetBackgroundColour(wx.NamedColor("WHITE"))

        self.center = outputfile.hierarchy.findMax("Density")[1]#[0.5,0.5,0.5]
        #self.center = [0.5, 0.5, 0.5]
        self.axis = axis

        self.makePlot()

        be.Initialize(canvas=FigureCanvas)
        self.figure_canvas = be.engineVals["canvas"](self, -1, self.figure)

        self.axes.set_xticks(())
        self.axes.set_yticks(())
        self.axes.set_ylabel("")
        self.axes.set_xlabel("")
        
        # Note that event is a MplEvent
        self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
        #self.figure_canvas.Bind(wx.EVT_ENTER_WINDOW, self.ChangeCursor)
        self.figure_canvas.Bind(wx.EVT_LEFT_UP, self.ClickProcess)
        self.figure_canvas.Bind(wx.EVT_CONTEXT_MENU, self.OnShowContextMenu)

        self.popupmenu = wx.Menu()
        self.cmapmenu = wx.Menu()
        self.popupmenu.AppendMenu(-1, "Color Map", self.cmapmenu)
        for cmap in be.matplotlib.cm.cmapnames:
            item = self.cmapmenu.AppendRadioItem(-1, cmap)
            self.Bind(wx.EVT_MENU, self.OnColorMapChoice, item)

        # Now we pre-select the current cmap
        ccmap_name = self.plot.colorbar.cmap.name
        self.cmapmenu.FindItemById(self.cmapmenu.FindItem(ccmap_name)).Check(True)

        self.widthSlider = wx.Slider(self, -1, wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        self.vals = na.logspace(log10(25*outputfile.hierarchy.getSmallestDx()),0,201)
        self.widthSlider.SetRange(0, 200)
        self.widthSlider.SetTickFreq(1,1)
        self.widthSlider.SetValue(200)

        self.widthBox = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.widthBox.SetValue("1.0")

        self.choices = outputfile.units.keys()
        self.choices.sort()

        self.unitList = wx.Choice(self, choices=self.choices)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.figure_canvas, 1, wx.EXPAND)

        self.ControlSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer.Add(self.ControlSizer, 0, wx.EXPAND)

        self.InputControl = wx.BoxSizer(wx.VERTICAL)
        self.InputControl.AddSpacer(10)
        self.InputControl.Add(self.widthBox, 0, 0, 0)
        self.InputControl.AddSpacer(10)
        self.InputControl.Add(self.unitList, 0, wx.EXPAND, 0)
        self.InputControl.AddSpacer(10)

        self.ControlSizer.AddSpacer(10)
        self.ControlSizer.Add(self.InputControl, 0, 0, 0)
        self.ControlSizer.AddSpacer(10)
        self.SliderSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SliderSizer.Add(self.widthSlider, 1,  wx.ALIGN_CENTER, 0)
        self.ControlSizer.Add(self.SliderSizer, 1, wx.GROW | wx.ALIGN_CENTER, 0)
        self.ControlSizer.AddSpacer(10)

        self.SetSizer(self.sizer)
        self.Fit()

        self.unitList.Bind(wx.EVT_CHOICE, self.UpdateUnit)
        self.widthSlider.Bind(wx.EVT_SCROLL, self.UpdateTextFromScrollEvent)
        self.widthSlider.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.UpdateWidth)
        self.widthBox.Bind(wx.EVT_TEXT_ENTER, self.UpdateWidthFromText)

    def OnShowContextMenu(self, event):
        pos = event.GetPosition()
        pos = self.figure_canvas.ScreenToClient(pos)
        self.figure_canvas.PopupMenu(self.popupmenu, pos)

    def OnColorMapChoice(self, event):
        item = self.cmapmenu.FindItemById(event.GetId())
        text = item.GetText()
        self.plot.set_cmap(text)
        self.UpdateWidth()

    def ClickProcess(self, event):
        xp, yp = event.X, event.Y
        dx = (self.plot.xlim[1] - self.plot.xlim[0])/self.plot.pix[0]
        dy = (self.plot.ylim[1] - self.plot.ylim[0])/self.plot.pix[1]
        if event.AltDown():
            l, b, width, height = self.figure.axes[0].bbox.get_bounds()
            #print l, b, width, height
            xp, yp = event.X, event.Y
            x = (dx * (xp-l)) + self.plot.xlim[0]
            y = (dy * (yp-b)) + self.plot.ylim[0]
            print "Setting",x,y
            self.center[yt.lagos.x_dict[self.axis]] = x
            self.center[yt.lagos.y_dict[self.axis]] = y
            self.UpdateWidth(self.widthSlider.GetValue())
        elif event.ShiftDown():
            if self.AmDrawingCircle:
                dc=wx.ClientDC(self.figure_canvas)
                b=dc.GetBrush()
                b.SetStyle(wx.TRANSPARENT)
                dc.SetBrush(b)
                rp = sqrt((xp-self.x1)**2.0 + \
                          (yp-self.y1)**2.0)
                r = sqrt(((xp-self.x1)*dx)**2.0 + \
                         ((yp-self.y1)*dy)**2.0)
                dc.DrawCircle(self.x1,self.y1,rp)
                unitname = self.choices[self.unitList.GetSelection()]
                unit = self.outputfile[unitname]
                print "R: %0.5e %s (%0.9e, %0.9e, %0.9e)" % (r*unit, unitname,
                    self.center[0], self.center[1], self.center[2])
                self.AmDrawingCircle = False
            else:
                self.x1 = xp
                self.y1 = yp
                self.AmDrawingCircle = True
            
    def UpdateUnit(self, event):
        pos = self.widthSlider.GetValue()
        self.UpdateTextFromScroll(pos)

    def UpdateTextFromScrollEvent(self, event):
        self.UpdateTextFromScroll(event.GetPosition())

    def UpdateTextFromScroll(self, pos):
        unitname = self.choices[self.unitList.GetSelection()]
        unit = self.outputfile[unitname]
        val = self.vals[pos] * unit
        #print unit, val
        self.widthBox.SetValue("%0.5e" % (val))

    def UpdateWidthFromText(self, event):
        val = float(event.GetString())
        if val < 0:
            return
        unitname = self.choices[self.unitList.GetSelection()]
        unit = self.outputfile[unitname]
        # Now we figure out the closest slider tick
        dx = (log10(self.vals[0])-log10(self.vals[-1]))/201
        self.widthSlider.SetValue(200-int(log10(val)/dx))
        #print dx, val, log10(val)/dx, int(log10(val)/dx)
        self.plot.set_width(val, unit)
        self.UpdateCanvas()

    def UpdateWidth(self, pos = None):
        val = float(self.widthBox.GetValue())
        unitname = self.choices[self.unitList.GetSelection()]
        unit = self.outputfile[unitname]
        #print unit, val, self.unitList.GetSelection()
        self.plot.set_width(val, unit)
        self.UpdateCanvas()

    def ChangeCursor(self, event):
        self.figure_canvas.SetCursor(wx.StockCursor(wx.CURSOR_BULLSEYE))

    def UpdateStatusBar(self, event):
        if event.inaxes:
            xp, yp = event.xdata, event.ydata
            dx = abs(self.plot.xlim[0] - self.plot.xlim[1])/self.plot.pix[0]
            dy = abs(self.plot.ylim[0] - self.plot.ylim[1])/self.plot.pix[1]
            x = (dx * xp) + self.plot.xlim[0]
            y = (dy * yp) + self.plot.ylim[0]
            self.statusBar.SetStatusText("x=%0.12e y=%0.12e z=%0.12e" %
                                           (x,y,self.GetDataValue(xp,yp)),
                                           0)

    def GetDataValue(self, x, y):
        return self.plot.image._A[int(y), int(x)]

    def fulldomain(self, *args):
        self.plot.set_width(1,"1")
        self.UpdateCanvas()

    def set_width(self, *args):
        w, u = Toolbars.ChooseWidth(self.outputfile)
        self.plot.set_width(w,u)
        self.UpdateCanvas()

    def set_zlim(self, *args):
        zmin, zmax = Toolbars.ChooseLimits(self.plot)
        self.plot.set_zlim(zmin,zmax)
        self.UpdateCanvas()

    def redraw(self, *args):
        self.UpdateCanvas()

    def save(self, *args):
        pass

    def switch_field(self, *args):
        field = Toolbars.ChooseField(self.outputfile)
        self.plot.switch_z(field)
        self.UpdateCanvas()

    def UpdateCanvas(self, *args):
        self.figure_canvas.draw()

class ReasonSlicePlotPanel(ReasonVMPlotPanel):
    def makePlot(self):
        self.data = self.outputfile.hierarchy.slice(self.axis, 
                                    self.center[self.axis], self.field, self.center)
        self.plot = be.SlicePlot(self.data, self.field, figure=self.figure, axes=self.axes)
        
class ReasonProjPlotPanel(ReasonVMPlotPanel):
    def makePlot(self):
        self.data = self.outputfile.hierarchy.proj(self.axis, 
                                   self.field, weightField=None, center=self.center)
        self.plot = be.ProjectionPlot(self.data, self.field, figure=self.figure, axes=self.axes)
        
class ReasonFidoSelector(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: ReasonFidoSelector.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.panel_4 = wx.Panel(self, -1)
        self.panel_5 = wx.Panel(self.panel_4, -1)
        self.FidoOutputs = wx.TreeCtrl(self.panel_4, -1, style=wx.TR_HAS_BUTTONS|wx.TR_DEFAULT_STYLE|wx.SUNKEN_BORDER)
        self.label_1 = wx.StaticText(self.panel_5, -1, "label_1")

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        pass
        # begin wxGlade: ReasonFidoSelector.__set_properties
        self.SetTitle("frame_4")
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: ReasonFidoSelector.__do_layout
        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_6.Add(self.FidoOutputs, 1, wx.EXPAND, 0)
        sizer_7.Add(self.label_1, 0, 0, 0)
        self.panel_5.SetSizer(sizer_7)
        sizer_6.Add(self.panel_5, 1, wx.EXPAND, 0)
        self.panel_4.SetSizer(sizer_6)
        sizer_5.Add(self.panel_4, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_5)
        sizer_5.Fit(self)
        self.Layout()

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

class ReasonApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        frame_1 = ReasonMainWindow(None, -1)
        self.SetTopWindow(frame_1)
        frame_1.Show()
        return True

if __name__ == "__main__":
    app = ReasonApp(0)
    app.MainLoop()
