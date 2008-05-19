"""
Notebook pages and main plotpanel classes.

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

_WelcomeMessage = \
"""
Welcome to Reason.

This GUI was designed as a front-end to a pre-existing toolkit called yt.  As such, it is feature-incomplete -- there's so much more that you can do with the command line.  That aside, there's some fun stuff you can do just in here.

To access the various plotting functions, right click on a dataset in the tree over to the left.  Then you can right-click in the window to shuffle around the center, to change the color map, to save it, and so on.  And then give a shot to shift-clicking in two places -- the first will define the center, the second will define the outer-edge.  Then you can phase plot the sphere!

If you experience any troubles, drop me a line (matt@yt.spacepope.org) or just fill out a ticket at yt.spacepope.org.
"""

class PlotPanel(wx.Panel):
    def __init__(self, *args, **kwds):
        # begin wxGlade: ReasonVMPlotFrame.__init__
        self.parent = kwds["parent"]
        self.CurrentlyResizing = False
        wx.Panel.__init__(self, *args, **kwds)
        self._link_plots = True

        self._AUI_NOTEBOOK = wx.NewId()
        self.nb = wx.aui.AuiNotebook(self, self._AUI_NOTEBOOK)
        welcomeMessage = wx.StaticText(self.nb, -1, _WelcomeMessage)
        self.nb.AddPage(welcomeMessage, "Welcome!")

        # Good thing the aui module breaks naming convention
        self.nb.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSE,
                     self.OnPageClosed)
        self.nb.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED,
                     self.OnPageChanged)
        self.__set_properties()
        self.__do_layout()

    def OnPageClosed(self, event):
        deletedPage = event.Selection
        Publisher().sendMessage("page_deleted",deletedPage)

    def OnPageChanged(self, event):
        page = self.nb.GetPage(event.Selection)
        self.UpdateSubscriptions(page)
        Publisher().sendMessage("page_changed",page)

    def UpdateSubscriptions(self, page):
        pairs = [("width","ChangeWidthFromMessage"),
                 ("field","ChangeFieldFromMessage"),
                 ("limits","ChangeLimitsFromMessage"),
                 ("center","ChangeCenterFromMessage"),
                 ("cmap","ChangeColorMapFromMessage")]
        for m,f in pairs:
            Publisher().unsubAll(("viewchange",m))
        if not hasattr(page,'CreationID'): return
        cti = page.CreationID
        toSubscribe = []
        if self._link_plots:
            for i in range(self.nb.GetPageCount()):
                p = self.nb.GetPage(i)
                try:
                    if p.CreationID == cti: toSubscribe.append(p)
                except AttributeError:
                    continue
        else: toSubscribe.append(page)
        for p in toSubscribe:
            Publisher().subscribe(self.MessageLogger)
            for m, f in pairs:
                ff = getattr(p,f)
                Publisher().subscribe(ff, ('viewchange',m))
        page.UpdateCanvas()

    def MessageLogger(self, message):
        print message

    def __set_properties(self):
        pass

    def __do_layout(self):
        PlotPanelSizer = wx.BoxSizer(wx.VERTICAL)
        PlotPanelSizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(PlotPanelSizer)
        self.Layout()

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

    # For the following, we need to figure out which windows to apply the
    # action to before we apply it.
    # This is probably best accomplished through pubsub.

    def AddPlot(self, page, title, ID):
        self.nb.AddPage(page, title)
        self.nb.SetSelection(self.nb.GetPageCount()+1)
        # Now we subscribe to events with this ID

    def GetCurrentOutput(self):
        return self.nb.GetPage(self.nb.Selection).outputfile

    def OnCallSwitchField(self, event):
        self.nb.GetPage(self.nb.Selection).switch_field(event)

    def OnCallSetWidth(self, event):
        self.nb.GetPage(self.nb.Selection).set_width(event)

    def OnCallRedraw(self, event):
        self.nb.GetPage(self.nb.Selection).UpdateCanvas()

    def OnCallZoomTop(self, event):
        self.nb.GetPage(self.nb.Selection).fulldomain(event)

    def OnCallSetZLim(self, event):
        self.nb.GetPage(self.nb.Selection).set_zlim(event)

    def OnCallViewPF(self, event):
        of = self.GetCurrentOutput()

class PlotPage(wx.Panel):
    plot = None
    def __init__(self, parent, status_bar, mw=None, CreationID = -1):
        wx.Panel.__init__(self, parent)

        self.CreationID = CreationID
        self.parent = parent
        self.mw = mw

        self.figure = be.matplotlib.figure.Figure((1,1))
        self.axes = self.figure.add_subplot(111)
        self.status_bar = status_bar

        self.SetupControls()
        self.SetupFigure()
        self.SetupMenu()
        self.DoLayout()

    def set_zlim(self, *args):
        zmin, zmax = Toolbars.ChooseLimits(self.plot)
        self.ChangeLimits(zmin, zmax)
        Publisher().sendMessage(("viewchange","limits"),(zmin,zmax))

    def redraw(self, *args):
        self.UpdateCanvas()

    # These should be done with better DRY style
    def OnEditTitle(self, event):
        curTitle = self.figure.axes[-1].title.get_text()
        dlg = wx.TextEntryDialog(
                self, 'New title?',
                'Change Title?', curTitle)
        if dlg.ShowModal() == wx.ID_OK:
            self.figure.axes[0].set_title(dlg.GetValue())
        dlg.Destroy()
        self.figure_canvas.draw()

    def OnEditLegend(self, event):
        try:
            curLegend = self.plot.colorbar.ax.yaxis.get_label().get_text()
        except:
            curLegend = ""
        dlg = wx.TextEntryDialog(
                self, 'New legend?',
                'Change Legend?', curLegend)
        if dlg.ShowModal() == wx.ID_OK:
            self.plot.colorbar.set_label(dlg.GetValue())
        dlg.Destroy()
        self.figure_canvas.draw()

    def OnEditXLabel(self, event):
        try:
            curXLabel = self.plot.axis.xaxis.get_label().get_text()
        except:
            curXLabel = ""
        dlg = wx.TextEntryDialog(
                self, 'New x-label?',
                'Change x-label?', curXLabel)
        if dlg.ShowModal() == wx.ID_OK:
            self.plot.axes.set_xlabel(dlg.GetValue())
        dlg.Destroy()
        self.figure_canvas.draw()

    def OnEditYLabel(self, event):
        try:
            curYLabel = self.plot.axis.yaxis.get_label().get_text()
        except:
            curYLabel = ""
        dlg = wx.TextEntryDialog(
                self, 'New y-label?',
                'Change y-label?', curYLabel)
        if dlg.ShowModal() == wx.ID_OK:
            self.plot.axes.set_ylabel(dlg.GetValue())
        dlg.Destroy()
        self.figure_canvas.draw()

    def OnShowContextMenu(self, event):
        pos = event.GetPosition()
        pos = self.figure_canvas.ScreenToClient(pos)
        self.ContextMenuPosition = pos
        self.figure_canvas.PopupMenu(self.popupmenu, pos)

    def OnColorMapChoice(self, event):
        item = self.cmapmenu.FindItemById(event.GetId())
        text = item.GetText()
        Publisher().sendMessage(("viewchange","cmap"),text)

    def ChangeLimits(self, zmin, zmax):
        self.plot.set_zlim(zmin,zmax)
        self.figure_canvas.draw()
        # We don't call update canvas

    def ChangeColorMapFromMessage(self, message):
        cmap = message.data
        self.plot.set_cmap(cmap)
        self.UpdateCanvas()

    def ChangeLimitsFromMessage(self, message):
        zmin, zmax = message.data
        self.ChangeLimits(zmin,zmax)
        self.UpdateCanvas()

    def ChangeFieldFromMessage(self, message):
        pass

    def ChangeCenterFromMessage(self, message):
        pass

    def ChangeWidthFromMessage(self, message):
        pass

    def SetupMenu(self):
        self.popupmenu = wx.Menu()
        self.cmapmenu = wx.Menu()
        self.editprops = wx.Menu()
        self.popupmenu.AppendMenu(-1, "Color Map", self.cmapmenu)
        cmapnames = sorted(color_maps)
        for cmap in cmapnames:
            item = self.cmapmenu.AppendRadioItem(-1, cmap)
            self.Bind(wx.EVT_MENU, self.OnColorMapChoice, item)
        self.popupmenu.AppendMenu(-1, "Edit Properties", self.editprops)
        dd = [('Edit Legend',self.OnEditLegend),
              ('Edit Title',self.OnEditTitle),
              ('Edit x-label',self.OnEditXLabel),
              ('Edit y-label',self.OnEditYLabel)]
        for title, func in dd:
            ii = self.editprops.Append(-1, title)
            self.Bind(wx.EVT_MENU, func, ii)
        # Now we pre-select the current cmap
        try:
            ccmap_name = self.plot.colorbar.cmap.name
        except AttributeError:
            ccmap_name = "jet"
        self.cmapmenu.FindItemById(self.cmapmenu.FindItem(ccmap_name)).Check(True)

    def makePlot(self):
        pass

    def UpdateStatusBar(self, event):
        pass

    def ClickProcess(self, event):
        pass

    def SaveImage(self):
        self.plot._generate_prefix('Plot')
        wildcard = "PNG (*png)|*.png"
        dlg = wx.FileDialog( \
            self, message="Save Image As ...", defaultDir=os.getcwd(), \
            defaultFile=self.plot.prefix, wildcard=wildcard, style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            orig_size = self.figure.get_size_inches()
            self.figure.set_size_inches((10,8))
            self.figure_canvas.print_figure(path,format='png')
            self.figure.set_size_inches(orig_size)
        dlg.Destroy()

    def SetupFigure(self):
        self.makePlot()
        be.Initialize(canvas=FigureCanvas)
        self.figure_canvas = be.engineVals["canvas"](self, -1, self.figure)
        self.figure_canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
        self.figure_canvas.Bind(wx.EVT_LEFT_DOWN, self.ClickProcess)
        self.figure_canvas.Bind(wx.EVT_CONTEXT_MENU, self.OnShowContextMenu)

class VMPlotPage(PlotPage):
    def __init__(self, parent, status_bar, outputfile, axis, field="Density",
                 weight_field = None, mw=None, CreationID = -1):
        self.outputfile = outputfile
        self.field = field
        self.weight_field = weight_field
        self.axis = axis
        self.center = [0.5, 0.5, 0.5]
        self.AmDrawingCircle = False
        self.circles = []

        if ytcfg.getboolean("reason","centeronmax"):
            self.center = outputfile.hierarchy.findMax("Density")[1]

        PlotPage.__init__(self, parent, status_bar, mw, CreationID)
        self.SetBackgroundColour(wx.NamedColor("WHITE"))
        self.UpdateWidth()

    def SetupMenu(self):
        PlotPage.SetupMenu(self)
        self.popupmenu.AppendSeparator()
        centerOnMax = self.popupmenu.Append(-1, "Center on max")
        centerHere = self.popupmenu.Append(-1, "Center here")
        self.popupmenu.AppendSeparator()
        gridBoundaries = self.popupmenu.AppendCheckItem(-1, "Show Grid Boundaries")
        gridBoundaries.Check(False)
        velocities = self.popupmenu.AppendCheckItem(-1, "Show Velocities")
        velocities.Check(False)
        self.popupmenu.AppendSeparator()
        self.take_log_menu = self.popupmenu.AppendCheckItem(-1, "Take Log")
        self.take_log_menu.Check(self.plot.log_field)
        fullDomain = self.popupmenu.Append(-1, "Zoom Top")

        self.Bind(wx.EVT_MENU, self.OnCenterOnMax, centerOnMax)
        self.Bind(wx.EVT_MENU, self.OnCenterHere, centerHere)
        self.Bind(wx.EVT_MENU, self.show_grid_boundaries, gridBoundaries)
        self.Bind(wx.EVT_MENU, self.show_velocities, velocities)
        self.Bind(wx.EVT_MENU, self.take_log, self.take_log_menu)
        self.Bind(wx.EVT_MENU, self.fulldomain, fullDomain)

    def SetupFigure(self):
        PlotPage.SetupFigure(self)

        self.axes.set_xticks(())
        self.axes.set_yticks(())
        self.axes.set_ylabel("")
        self.axes.set_xlabel("")

    def SetupControls(self):

        self.widthSlider = wx.Slider(self, -1, wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        self.vals = na.logspace(log10(25*self.outputfile.hierarchy.get_smallest_dx()),0,201)
        self.widthSlider.SetRange(0, 200)
        self.widthSlider.SetTickFreq(1,1)
        self.widthSlider.SetValue(200)

        self.widthBox = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.widthBox.SetValue("1.0")

        self.choices = self.outputfile.units.keys()
        self.choices.sort()

        self.unitList = wx.Choice(self, choices=self.choices)

        self.unitList.Bind(wx.EVT_CHOICE, self.UpdateUnit)
        self.widthSlider.Bind(wx.EVT_SCROLL, self.UpdateTextFromScrollEvent)
        self.widthSlider.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.UpdateWidth)
        self.widthBox.Bind(wx.EVT_TEXT_ENTER, self.UpdateWidthFromText)

        self.unitList.SetSelection(0)

    def DoLayout(self):

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

        self.ChangeWidth(1,"1")

    def ChangeColorMapFromMessage(self, message):
        PlotPage.ChangeColorMapFromMessage(self, message)
        self.UpdateWidth()

    def OnCenterOnMax(self, event):
        v, c = self.outputfile.h.findMax("Density")
        Publisher().sendMessage(('viewchange','center'), c)
        self.UpdateWidth()

    _grid_boundaries_cbid = None
    def show_grid_boundaries(self, event):
        if self._grid_boundaries_cbid is not None:
            self.plot.remove_callback(self._grid_boundaries_cbid)
            self._grid_boundaries_cbid = None
        else:
            self._grid_boundaries_cbid = \
                self.plot.add_callback(raven.be.gridBoundaryCallback())

    _velocities_cbid = None
    def show_velocities(self, event):
        if self._velocities_cbid is not None:
            self.plot.remove_callback(self._velocities_cbid)
            self._velocities_cbid = None
        else:
            xv = "%s-velocity" % (lagos.axis_names[lagos.x_dict[self.axis]])
            yv = "%s-velocity" % (lagos.axis_names[lagos.y_dict[self.axis]])
            self._velocities_cbid = \
                self.plot.add_callback(raven.be.quiverCallback(xv,yv,self.axis,20))

    def OnCenterHere(self, event):
        xp, yp = self.ContextMenuPosition
        x, y = self.ConvertPositionToDataPosition(xp, yp)
        print "CENTER HERE:", xp, yp, x, y
        if x == None or y == None: return
        self.ChangeCenter(x,y)
        self.UpdateWidth()

    def ChangeCenterFromMessage(self, message):
        #print "Hey, got message", message
        x, y, z = message.data
        # We are dealing with a pass-by-reference center
        self.center[0] = x
        self.center[1] = y
        self.center[2] = z
        self.data.reslice(self.center[self.axis])
        # If we do this, we get 3x as many updates as we want
        #self.UpdateWidth()

    def ConvertPositionToDataPosition(self, xp, yp):
        #print "CONVERT", xp, yp
        #if not self.figure.axes[0].in_axes(xp,yp): return None, None
        #xp, yp = self.figure.axes[0].transData.inverse_xy_tup((xp,yp))
        dx = (self.plot.xlim[1] - self.plot.xlim[0])/self.plot.pix[0]
        dy = (self.plot.ylim[1] - self.plot.ylim[0])/self.plot.pix[1]
        l, b, width, height = self.figure.axes[0].bbox.get_bounds()
        x = self.plot.xlim[0] + (dx * (xp-l))
        y = self.plot.ylim[1] - (dy * (yp-b))
        return x, y

    def ChangeCenter(self, x, y):
        newCenter = self.center[:]
        newCenter[lagos.x_dict[self.axis]] = x
        newCenter[lagos.y_dict[self.axis]] = y
        Publisher().sendMessage(('viewchange','center'), tuple(newCenter))
        self.UpdateWidth()

    def ClickProcess(self, event):
        xp, yp = event.X, event.Y
        if event.AltDown():
            x,y = self.ConvertPositionToDataPosition(xp, yp)
            self.ChangeCenter(x, y)
        elif event.ShiftDown():
            if self.AmDrawingCircle:
                dx = (self.plot.xlim[1] - self.plot.xlim[0])/self.plot.pix[0]
                dy = (self.plot.ylim[1] - self.plot.ylim[0])/self.plot.pix[1]
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
                cc = self.center[:]
                xd,yd = self.ConvertPositionToDataPosition(self.x1, self.y1)
                cc[lagos.x_dict[self.axis]] = xd
                cc[lagos.y_dict[self.axis]] = yd
                print "R: %0.5e %s (%0.9e, %0.9e, %0.9e)" % (r*unit, unitname,
                    cc[0], cc[1], cc[2])
                #print cc, r
                sphere = self.outputfile.hierarchy.sphere( \
                    cc, r, fields = ["Density"])
                self.mw._add_sphere("Sphere: %0.3e %s" \
                        % (r*unit, unitname), sphere)
                self.AmDrawingCircle = False
            else:
                self.x1 = xp
                self.y1 = yp
                self.AmDrawingCircle = True

    def ChangeWidthFromMessage(self, message):
        w,u=message.data
        # Now we need to update the text, units and bar
        wI = self.unitList.GetItems().index(u)
        self.unitList.SetSelection(wI)
        self.widthBox.SetValue("%0.5e" % (w))
        val = w/self.outputfile[u]
        dx = (log10(self.vals[0])-log10(self.vals[-1]))/201
        self.widthSlider.SetValue(200-int(log10(val)/dx))
        self.ChangeWidth(w,u)

    def ChangeWidth(self, width, unit):
        self.plot.set_width(width, unit)
        self.UpdateCanvas()

    def UpdateUnit(self, event=None):
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
        Publisher().sendMessage(('viewchange','width'), (val,unitname))

    def UpdateWidth(self, pos = None):
        val = float(self.widthBox.GetValue())
        unitname = self.choices[self.unitList.GetSelection()]
        unit = self.outputfile[unitname]
        Publisher().sendMessage(('viewchange','width'), (val,unitname))

    def ChangeCursor(self, event):
        self.figure_canvas.SetCursor(wx.StockCursor(wx.CURSOR_BULLSEYE))

    def UpdateStatusBar(self, event):
        #print event.x, event.y
        if event.inaxes:
            xp, yp = event.xdata, event.ydata
            dx = abs(self.plot.xlim[0] - self.plot.xlim[1])/self.plot.pix[0]
            dy = abs(self.plot.ylim[0] - self.plot.ylim[1])/self.plot.pix[1]
            x = (dx * xp) + self.plot.xlim[0]
            y = (dy * yp) + self.plot.ylim[0]
            xn = lagos.axis_names[lagos.x_dict[self.axis]]
            yn = lagos.axis_names[lagos.y_dict[self.axis]]
            self.status_bar.SetStatusText("%s = %0.12e" % (xn,x), 0)
            self.status_bar.SetStatusText("%s = %0.12e" % (yn,y), 1)
            self.status_bar.SetStatusText("v = %0.12e" % \
                                        (self.GetDataValue(xp,yp)), 2)

    def GetDataValue(self, x, y):
        return self.plot.image._A[int(y), int(x)]

    def take_log(self, event):
        self.plot.set_log_field(self.take_log_menu.IsChecked())
        self.UpdateCanvas()

    def fulldomain(self, *args):
        self.ChangeWidth(1,'1')
        Publisher().sendMessage(('viewchange','width'), (1,'1'))

    def set_width(self, *args):
        w, u = Toolbars.ChooseWidth(self.outputfile)
        self.ChangeWidth(w,u)
        Publisher().sendMessage(('viewchange','width'), (w,u))

    def ChangeFieldFromMessage(self, message):
        field = message.data
        self.ChangeField(field)

    def ChangeField(self, field):
        self.plot.switch_z(field)
        self.UpdateCanvas()
        self.take_log_menu.Check(self.plot.log_field)

    def switch_field(self, event):
        field = event.GetString()
        self.ChangeField(field)
        Publisher().sendMessage(('viewchange','field'), field)

    def UpdateCanvas(self, *args):
        if self.IsShown():
            self.plot._redraw_image()
            self.figure_canvas.draw()
        #else: print "Opting not to update canvas"

    def QueryFields(self):
        return QueryFields(self.outputfile)

class SlicePlotPage(VMPlotPage):
    def makePlot(self):
        self.data = self.outputfile.hierarchy.slice(self.axis,
                                    self.center[self.axis], self.field, self.center)
        self.plot = be.SlicePlot(self.data, self.field, figure=self.figure, axes=self.axes)

class ProjPlotPage(VMPlotPage):
    def makePlot(self):
        self.data = self.outputfile.hierarchy.proj(self.axis,
                                   self.field, weight_field=self.weight_field,
                                   center=self.center)
        self.plot = be.ProjectionPlot(self.data, self.field, figure=self.figure,
                                      axes=self.axes)

    def ChangeCenterFromMessage(self, message):
        x, y, z = message.data
        # We are dealing with a pass-by-reference center
        self.center[0] = x
        self.center[1] = y
        self.center[2] = z
        #self.UpdateWidth()

    def QueryFields(self):
        return [self.field]

class NewPhasePlotPage(PlotPage):
    def __init__(self, parent, status_bar, data_object, argdict, mw=None, CreationID = -1):
        self.data_object = data_object
        self.argdict = argdict

        PlotPage.__init__(self, parent, status_bar, mw, CreationID)
        self.UpdateAvailableFields()

    def SetupControls(self):
        self.ButtonPanel = wx.Panel(self, -1)
        self.sbox = wx.StaticBox(self.ButtonPanel, -1, "Field Selection")
        self.z_field = wx.Choice(self.ButtonPanel, -1, choices=[])
        self.add_field = wx.Button(self.ButtonPanel, -1, "Add Field")
        self.Bind(wx.EVT_BUTTON, self.AddField, self.add_field)
        self.Bind(wx.EVT_CHOICE, self.switch_z, self.z_field)

    def DoLayout(self):
        self.MainSizer = wx.BoxSizer(wx.VERTICAL)
        self.MainSizer.Add(self.figure_canvas, 1, wx.EXPAND)
        self.MainSizer.Add(self.ButtonPanel, 0, wx.EXPAND)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.z_field, 1, wx.EXPAND)
        sizer.Add(self.add_field, 0, wx.EXPAND)
        self.ButtonPanel.SetSizer(sizer)
        self.ButtonPanel.Layout()

        self.SetSizer(self.MainSizer)
        self.Layout()

    def AddField(self, event=None):
        pf2daf = Profile2DAddField(self.data_object, self)
        if pf2daf.ShowModal() == wx.ID_OK:
            argdict = pf2daf.return_argdict()
            if argdict is None:
                pf2daf.Destroy()
                return
            self.data.add_fields(**argdict)
            self.UpdateAvailableFields()
            self.z_field.SetSelection(
                    self.z_field.GetItems().index(argdict['fields']))
            self.switch_z()
        pf2daf.Destroy()

    def makePlot(self, event=None):
        self.data = lagos.BinnedProfile2D(self.data_object, **self.argdict)
        self.data.pf = self.data_object.pf
        self.data.hierarchy = self.data_object.pf.h
        # Now, unfortunately, we prompt for the field
        self.AddField(None)
        self.plot = be.NewPhasePlot(self.data, 
                                    [self.argdict['x_bin_field'],
                                     self.argdict['y_bin_field'],
                                     self.z_field.GetStringSelection()],
                                    figure=self.figure, axes=self.axes)

    def UpdateStatusBar(self, event):
        #print event.x, event.y
        if event.inaxes:
            if not hasattr(self.plot, 'pix'): return
            xp, yp = event.xdata, event.ydata
            xn = self.data.x_bin_field
            yn = self.data.y_bin_field
            vn = self.FieldZ.StringSelection
            self.status_bar.SetStatusText("%s = %0.5e" % (xn, xp), 0)
            self.status_bar.SetStatusText("%s = %0.5e" % (yn, yp), 1)
            self.status_bar.SetStatusText("%s = %0.5e" % \
                                        (vn, self.GetDataValue(xp,yp)), 2)

    def GetDataValue(self, x, y):
        xi = na.digitize(na.array([x]), self.plot.x_bins)-1
        yi = na.digitize(na.array([y]), self.plot.y_bins)-1
        return self.plot.vals[yi, xi]

    def UpdateAvailableFields(self):
        self.z_field.SetItems(sorted(self.QueryAvailableFields()))

    def QueryAvailableFields(self):
        fields = []
        for field in self.data.keys():
            if field != self.data.x_bin_field and \
               field != self.data.y_bin_field:
                fields.append(field)
        return fields

    def switch_z(self, event=None):
        if self.plot is None: return
        self.plot.switch_z(self.z_field.GetStringSelection())
        self.UpdateCanvas()

    def UpdateCanvas(self, *args):
        if self.IsShown():
            self.plot._redraw_image()
            self.figure_canvas.draw()
        #else: print "Opting not to update canvas"

    def ChangeLimits(self, zmin, zmax):
        print "Change Limits"
        self.plot.set_zlim(zmin,zmax)
        self.figure_canvas.draw()
        # We don't call update canvas

class PhasePlotPage(PlotPage):
    def __init__(self, parent, status_bar, data_object, mw=None, CreationID = -1):
        self.data_object = data_object

        PlotPage.__init__(self, parent, status_bar, mw, CreationID)

    def SetupControls(self):
        self.ButtonPanel = wx.Panel(self, -1)
        fs = self.QueryFields()
        self.FieldX = wx.Choice(self.ButtonPanel, -1, choices=fs, name="X")
        self.FieldX.SetSelection(0)
        self.FieldY = wx.Choice(self.ButtonPanel, -1, choices=fs, name="Y")
        self.FieldY.SetSelection(0)
        self.FieldZ = wx.Choice(self.ButtonPanel, -1, choices=fs, name="Z")
        self.FieldZ.SetSelection(0)
        self.FieldW = wx.Choice(self.ButtonPanel, -1, choices=fs, name="W")
        self.FieldW.SetSelection(0)
        self.FirePlot = wx.Button(self.ButtonPanel, label = "Make Plot")

        #self.Bind(wx.EVT_CHOICE, self.switch_x, self.FieldX)
        #self.Bind(wx.EVT_CHOICE, self.switch_y, self.FieldY)
        #self.Bind(wx.EVT_CHOICE, self.switch_z, self.FieldZ)
        #self.Bind(wx.EVT_CHOICE, self.switch_weight, self.FieldW)
        self.Bind(wx.EVT_BUTTON, self.fire_plot, self.FirePlot)

    def DoLayout(self):
        self.MainSizer = wx.BoxSizer(wx.VERTICAL)
        self.MainSizer.Add(self.figure_canvas, 1, wx.EXPAND)
        self.MainSizer.Add(self.ButtonPanel, 0, wx.EXPAND)

        self.ButtonSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.XSizer = wx.BoxSizer(wx.VERTICAL)
        self.YSizer = wx.BoxSizer(wx.VERTICAL)
        self.ZSizer = wx.BoxSizer(wx.VERTICAL)
        self.WSizer = wx.BoxSizer(wx.VERTICAL)


        self.ButtonSizer.AddSpacer(10)

        self.ButtonSizer.Add(self.XSizer, 1, wx.EXPAND)
        self.XSizer.Add(wx.StaticText(self.ButtonPanel, -1, "X Field"), 1, wx.EXPAND)
        self.XSizer.Add(self.FieldX, 1, wx.EXPAND)
        self.ButtonSizer.AddSpacer(20)

        self.ButtonSizer.Add(self.YSizer, 1, wx.EXPAND)
        self.YSizer.Add(wx.StaticText(self.ButtonPanel, -1, "Y Field"), 1, wx.EXPAND)
        self.YSizer.Add(self.FieldY, 1, wx.EXPAND)
        self.ButtonSizer.AddSpacer(20)

        self.ButtonSizer.Add(self.ZSizer, 1, wx.EXPAND)
        self.ZSizer.Add(wx.StaticText(self.ButtonPanel, -1, "Z Field"), 1, wx.EXPAND)
        self.ZSizer.Add(self.FieldZ, 1, wx.EXPAND)
        self.ButtonSizer.AddSpacer(20)

        self.ButtonSizer.Add(self.WSizer, 1, wx.EXPAND)
        self.WSizer.Add(wx.StaticText(self.ButtonPanel, -1, "Weight"), 1, wx.EXPAND)
        self.WSizer.Add(self.FieldW, 1, wx.EXPAND)

        self.ButtonSizer.AddSpacer(10)
        self.ButtonSizer.Add(self.FirePlot, 1, wx.EXPAND)
        self.ButtonSizer.AddSpacer(10)
        self.ButtonPanel.SetSizer(self.ButtonSizer)
        self.ButtonPanel.Layout()

        self.SetSizer(self.MainSizer)
        self.Layout()


    def QueryFields(self):
        return QueryFields(self.data_object)

    def makePlot(self, event=None):
        pass

    def fire_plot(self, event):
        if self.plot is not None:
            for ax in self.figure.axes[1:]: self.figure.delaxes(ax)
            self.axes.clear()
            self.figure.subplots_adjust()
            self.axes.set_xscale('linear')
            self.axes.set_xlim([0.0,1.0])
            self.axes.set_yscale('linear')
            self.axes.set_ylim([0.0,1.0])
            del self.plot
        X = self.FieldX.GetStringSelection()
        Y = self.FieldY.GetStringSelection()
        Z = self.FieldZ.GetStringSelection()
        W = self.FieldW.GetStringSelection()
        if len(W) == 0: W=None
        self.plot = be.PhasePlot(self.data_object, [X,Y,Z], weight = W,
                                 figure = self.figure, axes = self.axes)
        self.UpdateCanvas()

    def UpdateStatusBar(self, event):
        #print event.x, event.y
        if event.inaxes:
            if not hasattr(self.plot, 'pix'): return
            xp, yp = event.xdata, event.ydata
            xn = self.FieldX.StringSelection
            yn = self.FieldY.StringSelection
            vn = self.FieldZ.StringSelection
            self.status_bar.SetStatusText("%s = %0.5e" % (xn, xp), 0)
            self.status_bar.SetStatusText("%s = %0.5e" % (yn, yp), 1)
            self.status_bar.SetStatusText("%s = %0.5e" % \
                                        (vn, self.GetDataValue(xp,yp)), 2)

    def GetDataValue(self, x, y):
        xi = na.digitize(na.array([x]), self.plot.x_bins)-1
        yi = na.digitize(na.array([y]), self.plot.y_bins)-1
        return self.plot.vals[yi, xi]

    def switch_weight(self, event):
        if self.plot is None: return
        W = event.String
        if len(event.String) == 0: W = None
        self.plot.switch_weight(W)
        self.UpdateCanvas()

    def switch_x(self, event):
        if self.plot is None: return
        self.plot.switch_x(event.String)
        self.UpdateCanvas()

    def switch_y(self, event):
        if self.plot is None: return
        self.plot.switch_y(event.String)
        self.UpdateCanvas()

    def switch_z(self, event):
        if self.plot is None: return
        self.plot.switch_z(event.String)
        self.UpdateCanvas()

    def UpdateCanvas(self, *args):
        if self.IsShown():
            self.plot._redraw_image()
            self.figure_canvas.draw()
        #else: print "Opting not to update canvas"

    def ChangeLimitsFromMessage(self, message):
        zmin, zmax = message.data
        self.ChangeLimits(zmin,zmax)

    def ChangeLimits(self, zmin, zmax):
        print "Change Limits"
        self.plot.set_zlim(zmin,zmax)
        self.figure_canvas.draw()
        # We don't call update canvas

if __name__ == "__main__":
    app = ReasonApp(0)
    app.MainLoop()
