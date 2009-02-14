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

from yt.mods import *
pf = EnzoStaticOutput("/Users/matthewturk/Research/data/galaxy1200.dir/galaxy1200")

from enthought.traits.api import HasTraits, List, Instance, Str, Float, \
                    Any, Code, PythonValue
from enthought.traits.ui.api import Group, VGroup, HGroup, Tabbed, View, \
                    Item, ShellEditor, InstanceEditor, ListStrEditor, \
                    ListEditor, VSplit, VFlow, HSplit, VFold, ValueEditor, \
                    TreeEditor, TreeNode, RangeEditor
from enthought.traits.ui.menu import OKCancelButtons

class DataObjectList(HasTraits):
    data_objects = List(Str)

    view = View(
              Item('data_objects', show_label=False,
                   editor=ListStrEditor())
               )

    def _data_objects_default(self):
        return ['a','b','c']

class PlotSpec(HasTraits):
    width = Float(1.0)
    view = View(Item('width', editor=RangeEditor(),
                     show_label=False))

class PlotFrameTab(HasTraits):
    my_spec = Instance(PlotSpec)
    view = View(Item('my_spec', show_label=False, style='custom',
                     editor=InstanceEditor(editable=True)))

    def _my_spec_default(self):
        return PlotSpec()

class DataObject(HasTraits):
    name = Str

class ParameterFileNode(HasTraits):
    pf = Instance(EnzoStaticOutput)
    data_objects = List(Instance(DataObject))
    name = Str

    def _name_default(self):
        return str(self.pf)

class ParameterFileCollection(HasTraits):
    parameter_files = List(Instance(ParameterFileNode))
    name = Str

    def _parameter_files_default(self):
        gc = fido.GrabCollections()
        my_list = []
        for f in gc[0]:
            pf = EnzoStaticOutput(f)
            my_list.append(
                ParameterFileNode(pf=pf, 
                        data_objects = [DataObject(name='yo')]))
        return my_list

    def _name_default(self):
        gc = fido.GrabCollections()
        return str(gc[0])

class MainWindow(HasTraits):
    parameter_files = Instance(ParameterFileCollection)
    plot_frame_tabs = List(Instance(PlotFrameTab))
    shell = PythonValue

    def _shell_default(self):
        return globals()

    view = View(VSplit(
                    HSplit(
                       Item('parameter_files', 
                            height=700.0, width=100.0,
                            show_label=False,
                            editor = TreeEditor(editable=False,
                    nodes=[
                        TreeNode(node_for=[ParameterFileCollection],
                                 children='',
                                 label="=ParameterFiles"),
                        TreeNode(node_for=[ParameterFileCollection],
                                 children='parameter_files',
                                 label="name",
                                 view=View()),
                        TreeNode(node_for=[ParameterFileNode],
                                 children='data_objects',
                                 label="name",
                                 view=View()),
                        TreeNode(node_for=[DataObject],
                                 children='',
                                 label="name"),
                                ], show_icons=False),),
                       Item('plot_frame_tabs', style='custom',
                            editor=ListEditor(editor=InstanceEditor(editable=True),
                                              use_notebook=True),
                            show_label=False, height=700.0, width=800.0),
                    ),
                    HGroup(
                       Item('shell', editor=ShellEditor(share=True),
                            show_label=False, height=100.0),
                    ),
                ),
               resizable=True, width=0.9, height=0.9) 

    def my_select(self, ui):
        print "HI!"

    def _parameter_files_default(self):
        return ParameterFileCollection()

class YTScript(HasTraits):
    code = Code
    view = View(Item('code', show_label=False),
                height=0.8, width=0.8, resizable=True,
                buttons=OKCancelButtons)

class ObjectViewer(HasTraits):
    to_view=Any
    view = View(Item('to_view', editor=ValueEditor(), show_label=False),
                     resizable=True, height=0.8, width=0.8)

def view_object(obj):
    ObjectViewer(to_view=obj).edit_traits()

def run_script():
    my_script = YTScript()
    my_script.edit_traits()
    return my_script

dol = DataObjectList()
pft = [PlotFrameTab(), PlotFrameTab()]
ps = PlotSpec()
mw = MainWindow(data_object_list = dol,
                plot_frame_tabs = pft,
                plot_spec=ps)
mw.configure_traits()
