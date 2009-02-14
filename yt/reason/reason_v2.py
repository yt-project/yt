"""
New version of Reason, using a TraitsUI-based approach

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

from enthought.traits.api import HasTraits, List, Instance, Str, Float, \
                    Any
from enthought.traits.ui.api import Group, VGroup, HGroup, Tabbed, View, \
                    Item, ShellEditor, InstanceEditor, ListStrEditor, \
                    ListEditor, VSplit, VFlow, HSplit, VFold

class DataObjectList(HasTraits):
    data_objects = List(Str)

    view = View(
              Item('data_objects', show_label=False,
                   editor=ListStrEditor())
               )

    def _data_objects_default(self):
        return ['a','b','c']

class PlotFrameTab(HasTraits):
    my_name = Str('hi')
    view = View(Item('my_name', show_label=False))

class PlotSpec(HasTraits):
    width = Float(1.0)
    view = View(Item('width', show_label=False))

class MainWindow(HasTraits):
    data_object_list = Instance(DataObjectList)
    plot_frame_tabs = List(Instance(PlotFrameTab))
    plot_spec = Instance(PlotSpec)
    shell = Any

    view = View(VSplit(
                    HSplit(
                       Item('data_object_list', style='custom',
                            editor=InstanceEditor(editable=True),
                            show_label=False, height=700.0),
                       Item('plot_frame_tabs', style='custom',
                            editor=ListEditor(editor=InstanceEditor(editable=True),
                                              use_notebook=True),
                            show_label=False, height=700.0),
                       Item('plot_spec', style='custom',
                            editor=InstanceEditor(editable=True),
                            show_label=False, height=700.0),
                    ),
                    HGroup(
                       Item('shell', editor=ShellEditor(),
                            show_label=False, height=100.0),
                    ),
                ),
               resizable=True, width=0.9, height=0.9) 

dol = DataObjectList()
pft = [PlotFrameTab(), PlotFrameTab()]
ps = PlotSpec()
mw = MainWindow(data_object_list = dol,
                plot_frame_tabs = pft,
                plot_spec=ps)
mw.configure_traits()
