"""
Figure editors for the Traits GUI

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

import wx, sys, matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from enthought.traits.api import Any, Instance
from enthought.traits.ui.wx.editor import Editor
from enthought.traits.ui.wx.basic_editor_factory import BasicEditorFactory

from enthought.pyface.action.api import ActionController

from enthought.traits.ui.menu import \
    Menu, Action, Separator, OKCancelButtons, OKButton

from matplotlib.backend_bases import Event as MPLEvent

class _MPLFigureEditor(Editor):
    """ Snagged from Gael's tutorial """

    scrollable  = True
    mpl_control = Instance(FigureCanvas)

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.
        panel = wx.Panel(parent, -1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        self.mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(self.mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.value.canvas.SetMinSize((10,10))
        return panel

class MPLFigureEditor(BasicEditorFactory):
    klass = _MPLFigureEditor

class MPLAction(Action):
    event = Instance(MPLEvent)

class _MPLVMPlotEditor(_MPLFigureEditor, ActionController):

    def _create_canvas(self, parent):
        panel = _MPLFigureEditor._create_canvas(self, parent)
        self.mpl_control.mpl_connect("button_press_event", self.on_click)
        return panel

    def on_click(self, event):
        if not event.inaxes: return
        if event.button == 3:
            print "Clicky clicky"
            print event.xdata, event.ydata
            my_menu = Menu(MPLAction(name="Recenter", action="object.recenter",
                                     event=event),
                           MPLAction(name="Yo!", action="object.do_something",
                                     event=event))
            wxmenu = my_menu.create_menu(self.mpl_control, self)
            self.mpl_control.PopupMenuXY(wxmenu)

    def perform ( self, action ):
        """
        This is largely taken/modified from the TreeEditor _perform method.
        """
        object            = self.object
        method_name       = action.action
        info              = self.ui.info
        handler           = self.ui.handler
        event             = action.event

        if method_name.find( '.' ) >= 0:
            if method_name.find( '(' ) < 0:
                method_name += '(event)'
            try:
                eval( method_name, globals(),
                      { 'object':  object,
                        'editor':  self,
                        'info':    info,
                        'event':   event,
                        'handler': handler } )
            except:
                # fixme: Should the exception be logged somewhere?
                print sys.exc_info()
                
            return

        method = getattr( handler, method_name, None )
        if method is not None:
            method( info, object )
            return

        if action.on_perform is not None:
            action.on_perform( object )

class MPLVMPlotEditor(BasicEditorFactory):
    klass = _MPLVMPlotEditor

