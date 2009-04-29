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
#pf = EnzoStaticOutput("/Users/matthewturk/Research/data/galaxy1200.dir/galaxy1200")

from enthought.traits.api import \
    HasTraits, List, Instance, Str, Float, Any, Code, PythonValue, Int, CArray, \
    Property, Enum, cached_property, DelegatesTo, Callable, Array, \
    Button
from enthought.traits.ui.api import \
    Group, VGroup, HGroup, Tabbed, View, Item, ShellEditor, InstanceEditor, ListStrEditor, \
    ListEditor, VSplit, VFlow, HSplit, VFold, ValueEditor, TreeEditor, TreeNode, RangeEditor, \
    EnumEditor, Handler, Controller, DNDEditor
from enthought.traits.ui.menu import \
    Menu, Action, Separator, OKCancelButtons, OKButton
from enthought.pyface.action.api import \
    ActionController
from enthought.tvtk.pyface.scene_editor import SceneEditor
from enthought.tvtk.pyface.api import \
    DecoratedScene
from enthought.tvtk.pyface.scene_model import SceneModel
from enthought.traits.ui.wx.range_editor import SimpleSliderEditor

from plot_editors import Figure, MPLFigureEditor, MPLVMPlotEditor, Axes

from yt.raven.PlotTypes import VMPlot, ProjectionPlot, SlicePlot

import traceback
from tvtk_interface import \
    HierarchyImporter, YTScene

class PlotCreationHandler(Controller):
    main_window = Instance(HasTraits)
    pnode = Instance(HasTraits)

    format = Str
    plot_type = Any
    
    def close(self, info, is_ok):
        if not is_ok:
            super(Controller, self).close(info, True)
            return
        spt = self.plot_type(plot_spec=self.model, pf=self.pnode.pf,
                           name=self.format % (self.model.axis))
        self.pnode.data_objects.append(spt)
        self.main_window.plot_frame_tabs.append(spt)
        spt.plot

class VTKSceneCreationHandler(PlotCreationHandler):
    importer = Instance(HierarchyImporter)

    def close(self, info, is_ok):
        if is_ok: 
            yt_scene = YTScene(importer=self.importer,
                scene=SceneModel())
            spt = VTKDataObject(name = "VTK: %s" % self.pnode.pf,
                    scene=yt_scene.scene,
                    yt_scene=yt_scene)
            self.pnode.data_objects.append(spt)
            self.main_window.plot_frame_tabs.append(spt)
        super(Controller, self).close(info, True)
        return True


class DataObject(HasTraits):
    name = Str

class VTKDataObject(DataObject):
    yt_scene = Instance(YTScene)
    scene = DelegatesTo("yt_scene")
    add_contours = Button
    add_x_plane = Button
    add_y_plane = Button
    add_z_plane = Button
    edit_camera = Button
    edit_operators = Button
    edit_pipeline = Button
    center_on_max = Button
    operators = DelegatesTo("yt_scene")
    traits_view = View(
            Item("scene", editor = 
        SceneEditor(scene_class=DecoratedScene),
                    resizable=True, show_label=False),
            HGroup(Item("add_contours", show_label=False),
                   Item("add_x_plane", show_label=False),
                   Item("add_y_plane", show_label=False),
                   Item("add_z_plane", show_label=False),
                   Item("edit_camera", show_label=False),
                   Item("edit_operators", show_label=False),
                   Item("edit_pipeline", show_label=False),
                   Item("center_on_max", show_label=False),
                ),
            )

    operators_edit = View(
        Item("operators", style='custom', show_label=False,
             editor=ListEditor(editor=InstanceEditor(),
                               use_notebook=True),
              name="Edit Operators"),
        height=500.0, width=500.0, resizable=True)
    
    def _edit_camera_fired(self):
        self.yt_scene.camera_path.edit_traits()

    def _edit_operators_fired(self):
        self.edit_traits(view='operators_edit')

    def _edit_pipeline_fired(self):
        from enthought.tvtk.pipeline.browser import PipelineBrowser
        pb = PipelineBrowser(self.scene)
        pb.show()

    def _add_contours_fired(self):
        self.yt_scene.add_contour()

    def _add_x_plane_fired(self):
        self.yt_scene.add_x_plane()

    def _add_y_plane_fired(self):
        self.yt_scene.add_y_plane()

    def _add_z_plane_fired(self):
        self.yt_scene.add_z_plane()

    def _center_on_max_fired(self):
        self.yt_scene.do_center_on_max()

class ParameterFile(HasTraits):
    pf = Instance(EnzoStaticOutput)
    data_objects = List(Instance(DataObject))
    name = Str

    def _name_default(self):
        return str(self.pf)

    def do_slice(self):
        cons_view = View(
                Item('axis'), 
                Item('center'), 
                Item('field', editor=EnumEditor(name='field_list')),
                buttons=OKCancelButtons, title="Slicer: %s" % self.pf)
        ps = SlicePlotSpec(pf=self.pf)
        hand = PlotCreationHandler(main_window=mw, pnode=self, model=ps,
                                   plot_type=SlicePlotTab, format="Slice: %s")
        ps.edit_traits(cons_view, handler=hand)

    def do_proj(self):
        cons_view = View(
                Item('axis'), 
                Item('field', editor=EnumEditor(name='field_list')),
                Item('weight_field', editor=EnumEditor(name='none_field_list')),
                buttons=OKCancelButtons, title="Projector: %s" % self.pf)
        ps = ProjPlotSpec(pf=self.pf)
        hand = PlotCreationHandler(main_window=mw, pnode=self, model=ps,
                                   plot_type=ProjPlotTab, format="Proj: %s")
        ps.edit_traits(cons_view, handler=hand)

    def do_vtk(self):
        from tvtk_interface import HierarchyImporter, \
            HierarchyImportHandler
        importer = HierarchyImporter(pf=self.pf, max_level=self.pf.h.max_level)
        importer.edit_traits(handler = VTKSceneCreationHandler(
            main_window=mw, pnode=self, importer = importer))

class ParameterFileCollection(HasTraits):
    parameter_files = List(Instance(ParameterFile))
    name = Str
    collection = Any

    def _parameter_files_default(self):
        my_list = []
        for f in self.collection:
            try:
                pf = EnzoStaticOutput(f)
                my_list.append(
                    ParameterFile(pf=pf, 
                            data_objects = []))
            except IOError: pass
        return my_list

    def _name_default(self):
        return str(self.collection)

class ParameterFileCollectionList(HasTraits):
    parameter_file_collections = List(Instance(ParameterFileCollection))

    def _parameter_file_collections_default(self):
        return [ParameterFileCollection(collection=c)
                for c in fido.GrabCollections()]

class DataObjectList(HasTraits):
    data_objects = List(Str)

    traits_view = View(
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

class VMPlotSpec(HasTraits):
    pf = Instance(EnzoStaticOutput)
    field = Str('Density')
    field_list = Property(depends_on = 'pf')

    center = Array(shape=(3,), dtype='float64')
    axis = Enum(0,1,2)

    @cached_property
    def _get_field_list(self):
        fl = self.pf.h.field_list
        df = self.pf.h.derived_field_list
        fl.sort(); df.sort()
        return fl + df

    def _center_default(self):
        return self.pf.h.find_max("Density")[1]

class SlicePlotSpec(VMPlotSpec):
    pass

class ProjPlotSpec(VMPlotSpec):
    weight_field = Str("None")
    none_field_list = Property(depends_on = 'field_list')

    @cached_property
    def _get_none_field_list(self):
        return ["None"] + self.field_list

class VMPlotTab(PlotFrameTab):
    pf = Instance(EnzoStaticOutput)
    figure = Instance(Figure, args=())
    field = DelegatesTo('plot_spec')
    field_list = DelegatesTo('plot_spec')
    plot = Instance(VMPlot)
    axes = Instance(Axes)
    disp_width = Float(1.0)
    unit = Str('unitary')
    min_width = Property(Float, depends_on=['pf','unit'])
    max_width = Property(Float, depends_on=['pf','unit'])
    unit_list = Property(depends_on = 'pf')
    smallest_dx = Property(depends_on = 'pf')

    traits_view = View(VGroup(
            HGroup(Item('figure', editor=MPLVMPlotEditor(),
                     show_label=False)),
            HGroup(Item('disp_width',
                     editor=RangeEditor(format="%0.2e",
                        low_name='min_width', high_name='max_width',
                        mode='logslider', enter_set=True),
                     show_label=False, width=400.0),
                   Item('unit',
                      editor=EnumEditor(name='unit_list')),),
            HGroup(Item('field',
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
        return 50.0*self.smallest_dx*self.pf[self.unit]

    @cached_property
    def _get_max_width(self):
        return self.pf['unitary']*self.pf[self.unit]

    @cached_property
    def _get_smallest_dx(self):
        return self.pf.h.get_smallest_dx()

    @cached_property
    def _get_unit_list(self):
        return self.pf.units.keys()

    def _unit_changed(self, old, new):
        self.disp_width = self.disp_width * self.pf[new]/self.pf[old]

    def _disp_width_changed(self, old, new):
        self.plot.set_width(new, self.unit)
        self._redraw()

    def _redraw(self):
        self.figure.canvas.draw()

    def recenter(self, event):
        xp, yp = event.xdata, event.ydata
        dx = abs(self.plot.xlim[0] - self.plot.xlim[1])/self.plot.pix[0]
        dy = abs(self.plot.ylim[0] - self.plot.ylim[1])/self.plot.pix[1]
        x = (dx * xp) + self.plot.xlim[0]
        y = (dy * yp) + self.plot.ylim[0]
        xi = lagos.x_dict[self.axis]
        yi = lagos.y_dict[self.axis]
        cc = self.center[:]
        cc[xi] = x; cc[yi] = y
        self.plot.data.center = cc[:]
        self.plot.data.set_field_parameter('center', cc.copy())
        self.center = cc

class SlicePlotTab(VMPlotTab):
    plot_spec = Instance(SlicePlotSpec)

    axis = DelegatesTo('plot_spec')
    center = DelegatesTo('plot_spec')
    
    plot = Instance(SlicePlot)

    def _plot_default(self):
        coord = self.center[self.axis]
        sl = self.pf.h.slice(self.axis, coord, center=self.center[:])
        sp = SlicePlot(sl, self.field, self.figure, self.axes)
        self.figure.canvas.draw()
        return sp

    def _center_changed(self, old, new):
        #traceback.print_stack()
        if na.all(na.abs(old - new) == 0.0): return
        print na.abs(old-new)
        print "Re-slicing", old, new
        pp = self.center
        self.plot.data.reslice(pp[self.axis])
        self.plot._refresh_display_width()
        self.figure.canvas.draw()

class ProjPlotTab(VMPlotTab):
    plot_spec = Instance(ProjPlotSpec)

    axis = DelegatesTo('plot_spec')
    center = DelegatesTo('plot_spec')
    weight_field = DelegatesTo('plot_spec')

    plot = Instance(ProjectionPlot)

    def _plot_default(self):
        self.field = self.field[:]
        self.weight_field = self.weight_field[:]
        wf = self.weight_field
        if str(wf) == "None": wf = None
        proj = self.pf.h.proj(self.axis, self.field, wf,
                        center=self.center[:])
        pp = ProjectionPlot(proj, self.field, self.figure, self.axes)
        self.figure.canvas.draw()
        return pp

    def _center_changed(self, old, new):
        self.plot._refresh_display_width()

class SphereWrapper(DataObject):
    radius = Float
    unit = Str

class MainWindow(HasTraits):
    parameter_file_collections = Instance(ParameterFileCollectionList)
    parameter_files = Instance(ParameterFileCollection)
    plot_frame_tabs = List(Instance(DataObject))
    open_parameterfile = Button
    shell = PythonValue

    def _shell_default(self):
        return globals()
    notebook_editor = ListEditor(editor=InstanceEditor(editable=True),
                                 use_notebook=True)

    traits_view = View(VSplit(
                    HSplit(VGroup(
                       Item('parameter_file_collections', 
                            width=120.0, height=500.0,
                            show_label=False,
                            editor = TreeEditor(editable=False,
                    nodes=[
                        TreeNode(node_for=[ParameterFileCollectionList],
                                 children='parameter_file_collections',
                                 label="=Data Collections"),
                        TreeNode(node_for=[ParameterFileCollection],
                                 children='parameter_files',
                                 label="name",
                                 view=View()),
                        TreeNode(node_for=[ParameterFile],
                                 children='data_objects',
                                 label="name",
                                 menu = Menu(Action(name='Slice',
                                                    action='object.do_slice'),
                                             Action(name='Project',
                                                    action='object.do_proj'),
                                             Action(name='VTK',
                                                    action='object.do_vtk')),
                                 view=View()),
                        TreeNode(node_for=[DataObject],
                                 children='',
                                 label="name"),
                                ], show_icons=False),),
                        Item('open_parameterfile', show_label=False)),
                       Item('plot_frame_tabs', style='custom',
                            editor = notebook_editor,
                            show_label=False, height=500.0, width=500.0),
                    ),
                    HGroup(
                       Item('shell', editor=ShellEditor(share=True),
                            show_label=False, height=120.0),
                    ),
                ),
               resizable=True, width=800.0, height=660.0,
               title="reason v2 [prototype]")

    def _open_parameterfile_fired(self):
        print "OPENING"

    def _parameter_file_collections_default(self):
        return ParameterFileCollectionList()

class YTScript(HasTraits):
    code = Code
    traits_view = View(Item('code', show_label=False),
                       height=0.8, width=0.8, resizable=True,
                       buttons=OKCancelButtons)

class ObjectViewer(HasTraits):
    to_view=Any
    traits_view = View(
            Item('to_view', editor=ValueEditor(), show_label=False),
                     resizable=True, height=0.8, width=0.8)

def view_object(obj):
    ObjectViewer(to_view=obj).edit_traits()

def run_script():
    my_script = YTScript()
    my_script.edit_traits()
    return my_script

class event_mock(object):
    inaxes = True
    button = 3

dol = DataObjectList()
mw = MainWindow(plot_frame_tabs = [])
mw.configure_traits()

