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

from enthought.traits.api import \
    HasTraits, List, Instance, Str, Float, Any, Code, PythonValue, Int, CArray, \
    Property, Enum, cached_property
from enthought.traits.ui.api import \
    Group, VGroup, HGroup, Tabbed, View, Item, ShellEditor, InstanceEditor, ListStrEditor, \
    ListEditor, VSplit, VFlow, HSplit, VFold, ValueEditor, TreeEditor, TreeNode, RangeEditor, \
    EnumEditor
from enthought.traits.ui.menu import \
    Menu, Action, Separator
from enthought.traits.ui.menu import \
    OKCancelButtons

from enthought.traits.ui.wx.range_editor import SimpleSliderEditor

from plot_editors import Figure, MPLFigureEditor, Axes

from yt.raven.PlotTypes import VMPlot, ProjectionPlot, SlicePlot

class DataObject(HasTraits):
    name = Str

class ParameterFile(HasTraits):
    pf = Instance(EnzoStaticOutput)
    data_objects = List(Instance(DataObject))
    name = Str

    def _name_default(self):
        return str(self.pf)

    def do_slice(self):
        spt = SlicePlotTab(pf=self.pf, name="MySlice")
        self.data_objects.append(spt)
        mw.plot_frame_tabs.append(spt)
        spt.plot

class ParameterFileCollection(HasTraits):
    parameter_files = List(Instance(ParameterFile))
    name = Str

    def _parameter_files_default(self):
        gc = fido.GrabCollections()
        my_list = []
        for f in gc[0]:
            pf = EnzoStaticOutput(f)
            my_list.append(
                ParameterFile(pf=pf, 
                        data_objects = []))
        return my_list

    def _name_default(self):
        gc = fido.GrabCollections()
        return str(gc[0])

class DataObjectList(HasTraits):
    data_objects = List(Str)

    view = View(
              Item('data_objects', show_label=False,
                   editor=ListStrEditor())
               )

    def _data_objects_default(self):
        return ['a','b','c']

class PlotFrameTab(DataObject):
    figure = Instance(Figure)

class _LogFormat(str):
    def _convert_floats(self, other):
        args = []
        if not isinstance(other, types.TupleType):
            other = (other,)
        for arg in other:
            if isinstance(arg, types.FloatType):
                args.append(10**arg)
            else:
                args.append(arg)
        return tuple(args)

    def __mod__(self, other):
        args = self._convert_floats(other)
        return str.__mod__(self, tuple(args))

    def __rmod__(self, other):
        args = self._convert_floats(other)
        return str.__rmod__(self, tuple(args))

lf = _LogFormat("%0.2e")

class VMPlotTab(PlotFrameTab):
    
    pf = Instance(EnzoStaticOutput)
    figure = Instance(Figure, args=())
    plot = Instance(VMPlot)
    axes = Instance(Axes)
    disp_width = Float(1.0)
    unit = Str('unitary')
    field = Str('Density')
    min_width = Property(Float, depends_on=['pf','unit'])
    max_width = Property(Float, depends_on=['pf','unit'])
    unit_list = Property(depends_on = 'pf')
    field_list = Property(depends_on = 'pf')
    smallest_dx = Property(depends_on = 'pf')

    view = View(VGroup(
            HGroup(Item('figure', editor=MPLFigureEditor(),
                     show_label=False)),
            HGroup(Item('disp_width',
                     editor=RangeEditor(format="%0.2e",
                        low_name='min_width', high_name='max_width'),
                     show_label=False),
                   Item('unit',
                      editor=EnumEditor(name='unit_list')),
                   Item('field',
                      editor=EnumEditor(name='field_list')),
                )),
             resizable=True)

    def __init__(self, **traits):
        super(VMPlotTab, self).__init__(**traits)
        self.axes = self.figure.add_subplot(111, aspect='equal')

    def _field_changed(self, old, new):
        self.plot.switch_z(new)
        self._redraw()

    @cached_property
    def _get_min_width(self):
        print "HI!"
        return 50.0*self.smallest_dx*self.pf[self.unit]

    @cached_property
    def _get_max_width(self):
        print "OH NO"
        return self.pf['unitary']*self.pf[self.unit]

    @cached_property
    def _get_smallest_dx(self):
        print "SMALLEST"
        return self.pf.h.get_smallest_dx()

    @cached_property
    def _get_unit_list(self):
        return self.pf.units.keys()

    @cached_property
    def _get_field_list(self):
        fl = self.pf.h.field_list
        df = self.pf.h.derived_field_list
        fl.sort(); df.sort()
        return fl + df

    def _unit_changed(self, old, new):
        self.disp_width = self.disp_width * self.pf[new]/self.pf[old]

    def _disp_width_changed(self, old, new):
        self.plot.set_width(new, self.unit)
        self._redraw()

    def _redraw(self):
        self.figure.canvas.draw()

class SlicePlotTab(VMPlotTab):
    plot = Instance(SlicePlot)
    axis = Int(0)
    center = CArray(shape=(3,), dtype='float64')

    def _plot_default(self):
        coord = self.center[self.axis]
        sl = self.pf.h.slice(self.axis, coord, center=self.center)
        sp = SlicePlot(sl, self.field, self.figure, self.axes)
        self.figure.canvas.draw()
        return sp

    def _center_default(self):
        return self.pf.h.find_max("Density")[1]

class SphereWrapper(DataObject):
    radius = Float
    unit = Str

class MainWindow(HasTraits):
    parameter_files = Instance(ParameterFileCollection)
    plot_frame_tabs = List(Instance(PlotFrameTab))
    shell = PythonValue

    def _shell_default(self):
        return globals()

    view = View(VSplit(
                    HSplit(
                       Item('parameter_files', 
                            width=120.0, height=500.0,
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
                        TreeNode(node_for=[ParameterFile],
                                 children='data_objects',
                                 label="name",
                                 menu = Menu(Action(name='slice',
                                                    action='object.do_slice')),
                                 view=View()),
                        TreeNode(node_for=[DataObject],
                                 children='',
                                 label="name"),
                                ], show_icons=False),),
                       Item('plot_frame_tabs', style='custom',
                            editor=ListEditor(editor=InstanceEditor(editable=True),
                                              use_notebook=True),
                            show_label=False, height=500.0, width=500.0),
                    ),
                    HGroup(
                       Item('shell', editor=ShellEditor(share=True),
                            show_label=False, height=120.0),
                    ),
                ),
               resizable=True, width=800.0, height=660.0) 

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
pft = [SlicePlotTab(pf=pf)]
mw = MainWindow(data_object_list = dol,
                plot_frame_tabs = [])
mw.configure_traits()
