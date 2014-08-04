"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# Major library imports
from vm_panner import VariableMeshPanner
from numpy import linspace, meshgrid, pi, sin, mgrid, zeros

# Enthought library imports
from enthought.enable.api import Component, ComponentEditor, Window
from enthought.traits.api import HasTraits, Instance, Button, Any, Callable, \
        on_trait_change, Bool, DelegatesTo, List, Enum, Int, Property, Str, \
        cached_property
from enthought.traits.ui.api import \
        Group, View, Item, VGroup, InstanceEditor

# Chaco imports
from enthought.chaco.api import ArrayPlotData, jet, Plot, HPlotContainer, \
        ColorBar, DataRange1D, DataRange2D, LinearMapper, ImageData, \
        CMapImagePlot, OverlayPlotContainer
from enthought.chaco.tools.api import PanTool, ZoomTool, RangeSelection, \
        RangeSelectionOverlay, RangeSelection
from enthought.chaco.tools.image_inspector_tool import ImageInspectorTool, \
     ImageInspectorOverlay

if not hasattr(DataRange2D, "_subranges_updated"):
    print ("You'll need to add _subranges updated to enthought/chaco/data_range_2d.py")
    print ('Add this at the correct indentation level:')
    print ()
    print ('    @on_trait_change("_xrange.updated,_yrange.updated")')
    print ('    def _subranges_updated(self):')
    print ('        self.updated = True')
    print ()
    raise RuntimeError

# We like the algae colormap; for now we re-implement it here.

from enthought.chaco.api import \
    ColorMapper, \
    color_map_functions, color_map_dict, color_map_name_dict

def algae(range, **traits):
    _data = {'red':   ((0.0, 80/256., 80/256.),
                       (0.2, 0.0, 0.0),
                       (0.4, 0.0, 0.0),
                       (0.6, 256/256., 256/256.),
                       (0.95, 256/256., 256/256.),
                       (1.0, 150/256., 150/256.)),
             'green': ((0.0, 0/256., 0/256.),
                       (0.2, 0/256., 0/256.),
                       (0.4, 130/256., 130/256.),
                       (0.6, 256/256., 256/256.),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 80/256., 80/256.),
                       (0.2, 220/256., 220/256.),
                       (0.4, 0.0, 0.0),
                       (0.6, 20/256., 20/256.),
                       (1.0, 0.0, 0.0))}
    return ColorMapper.from_segment_map(_data, range=range, **traits)
color_map_functions.append(algae)
color_map_dict[algae] = "algae"
color_map_name_dict["algae"] = algae

class FunctionImageData(ImageData):
    # The function to call with the low and high values of the range.
    # It should return an array of values.
    func = Callable

    # A reference to a datarange
    data_range = Instance(DataRange2D)

    def __init__(self, **kw):
        # Explicitly call the AbstractDataSource constructor because
        # the ArrayDataSource ctor wants a data array
        ImageData.__init__(self, **kw)
        self.recalculate()

    @on_trait_change('data_range.updated')
    def recalculate(self):
        if self.func is not None and self.data_range is not None:
            newarray = self.func(self.data_range.low, self.data_range.high)
            ImageData.set_data(self, newarray)
        else:
            self._data = zeros((512,512),dtype=float)

    def set_data(self, *args, **kw):
        raise RuntimeError("Cannot set numerical data on a FunctionDataSource")

    def set_mask(self, mask):
        raise NotImplementedError

    def remove_mask(self):
        raise NotImplementedError

class ImagePixelizerHelper(object):
    index = None
    def __init__(self, panner, run_callbacks = False):
        self.panner = panner
        self.run_callbacks = run_callbacks

    def __call__(self, low, high):
        b = self.panner.set_low_high(low, high)
        if self.run_callbacks:
            self.panner._run_callbacks()
        if self.index is not None:
            num_x_ticks = b.shape[0] + 1
            num_y_ticks = b.shape[1] + 1
            xs = mgrid[low[0]:high[0]:num_x_ticks*1j]
            ys = mgrid[low[1]:high[1]:num_y_ticks*1j]
            self.index.set_data( xs, ys )
        return b

class ZoomedPlotUpdater(object):
    fid = None
    def __init__(self, panner, zoom_factor=4):
        """
        Supply this an a viewport_callback argument to a panner if you want to
        update a second panner in a smaller portion at higher resolution.  If
        you then set the *fid* property, you can also have it update a
        FunctionImageData datarange.  *panner* is the panner to update (not the
        one this is a callback to) and *zoom_factor* is how much to zoom in by.
        """
        self.panner = panner
        self.zoom_factor = zoom_factor

    def __call__(self, xlim, ylim):
        self.panner.xlim = xlim
        self.panner.ylim = ylim
        self.panner.zoom(self.zoom_factor)
        nxlim = self.panner.xlim
        nylim = self.panner.ylim
        if self.fid is not None:
            self.fid.data_range.set_bounds(
                (nxlim[0], nylim[0]), (nxlim[1], nylim[1]))

class VMImagePlot(HasTraits):
    plot = Instance(Plot)
    fid = Instance(FunctionImageData)
    img_plot = Instance(CMapImagePlot)
    panner = Instance(VariableMeshPanner)
    helper = Instance(ImagePixelizerHelper)
    fields = List

    def __init__(self, *args, **kwargs):
        super(VMImagePlot, self).__init__(**kwargs)
        self.add_trait("field", Enum(*self.fields))
        self.field = self.panner.field

    def _plot_default(self):
        pd = ArrayPlotData()
        plot = Plot(pd, padding = 0)
        self.fid._data = self.panner.buffer

        pd.set_data("imagedata", self.fid)

        img_plot = plot.img_plot("imagedata", colormap=algae,
                                 interpolation='nearest',
                                 xbounds=(0.0, 1.0),
                                 ybounds=(0.0, 1.0))[0]
        self.fid.data_range = plot.range2d
        self.helper.index = img_plot.index
        self.img_plot = img_plot
        return plot

    def _field_changed(self, old, new):
        self.panner.field = new
        self.fid.recalculate()

    def _fid_default(self):
        return FunctionImageData(func = self.helper)

    def _helper_default(self):
        return ImagePixelizerHelper(self.panner)

    def _panner_changed(self, old, new):
        index = self.helper.index
        self.helper = ImagePixelizerHelper(new)
        self.helper.index = index
        self.fid.func = self.helper
        self.fid.recalculate()

    def _fields_default(self):
        keys = []
        for field in self.panner.source.keys():
            if field not in ['px','py','pdx','pdy',
                             'pz','pdz','weight_field']:
                keys.append(field)
        return keys

class VariableMeshPannerView(HasTraits):

    plot = Instance(Plot)
    spawn_zoom = Button
    vm_plot = Instance(VMImagePlot)
    use_tools = Bool(True)
    full_container = Instance(HPlotContainer)
    container = Instance(OverlayPlotContainer)
    
    traits_view = View(
                    Group(
                        Item('full_container',
                             editor=ComponentEditor(size=(512,512)), 
                             show_label=False),
                        Item('field', show_label=False),
                        orientation = "vertical"),
                    width = 800, height=800,
                    resizable=True, title="Pan and Scan",
                    )

    def _vm_plot_default(self):
        return VMImagePlot(panner=self.panner)
    
    def __init__(self, **kwargs):
        super(VariableMeshPannerView, self).__init__(**kwargs)
        # Create the plot
        self.add_trait("field", DelegatesTo("vm_plot"))

        plot = self.vm_plot.plot
        img_plot = self.vm_plot.img_plot

        if self.use_tools:
            plot.tools.append(PanTool(img_plot))
            zoom = ZoomTool(component=img_plot, tool_mode="box", always_on=False)
            plot.overlays.append(zoom)
            imgtool = ImageInspectorTool(img_plot)
            img_plot.tools.append(imgtool)
            overlay = ImageInspectorOverlay(component=img_plot, image_inspector=imgtool,
                                            bgcolor="white", border_visible=True)
            img_plot.overlays.append(overlay)


        image_value_range = DataRange1D(self.vm_plot.fid)
        cbar_index_mapper = LinearMapper(range=image_value_range)
        self.colorbar = ColorBar(index_mapper=cbar_index_mapper,
                                 plot=img_plot,
                                 padding_right=40,
                                 resizable='v',
                                 width=30)

        self.colorbar.tools.append(
            PanTool(self.colorbar, constrain_direction="y", constrain=True))
        zoom_overlay = ZoomTool(self.colorbar, axis="index", tool_mode="range",
                                always_on=True, drag_button="right")
        self.colorbar.overlays.append(zoom_overlay)

        # create a range selection for the colorbar
        range_selection = RangeSelection(component=self.colorbar)
        self.colorbar.tools.append(range_selection)
        self.colorbar.overlays.append(
                RangeSelectionOverlay(component=self.colorbar,
                                      border_color="white",
                                      alpha=0.8, fill_color="lightgray"))

        # we also want to the range selection to inform the cmap plot of
        # the selection, so set that up as well
        range_selection.listeners.append(img_plot)

        self.full_container = HPlotContainer(padding=30)
        self.container = OverlayPlotContainer(padding=0)
        self.full_container.add(self.colorbar)
        self.full_container.add(self.container)
        self.container.add(self.vm_plot.plot)

class OutputSelector(HasTraits):
    outputs = List
    main_plot = Instance(VariableMeshPannerView)
    main_panner = Property(depends_on = "ds")
    weight_field = Str("Density")

    ds = Any
    source = Any
    axis = Int(0)

    traits_view = View(VGroup(
                        Item('output'),
                        Item('main_plot'),
                        )
                      )

    def __init__(self, **kwargs):
        super(OutputSelector, self).__init__(**kwargs)
        self.add_trait("output", Enum(*self.outputs))
        self.output = self.outputs[-1]
        self.main_plot

    def _output_default(self):
        return self.outputs[0]

    def _output_changed(self, old, new):
        # We get a string here
        import yt.mods
        self.ds = yt.mods.load(new, dataset_type="enzo_packed_3d")
        self.source = yt.mods.projload(self.ds, self.axis, "Density")
        self.main_panner.field = self.main_plot.vm_plot.field
        self.main_plot.panner = self.main_plot.vm_plot.panner = \
            self.main_plot.vm_plot.helper.panner = self.main_panner
        self.main_plot.vm_plot.field = self.main_panner.field

    def _main_plot_default(self):
        vmpv = VariableMeshPannerView(panner = self.main_panner)
        vmpv.vm_plot.helper.run_callbacks = True
        return vmpv

    @cached_property
    def _get_main_panner(self):
        return self.ds.image_panner(self.source, (512, 512), "Density")

def pan_and_scan_directory(dir_name):
    import glob, os
    fns = [ fn[:-10] for fn in
            glob.glob(os.path.join(dir_name, "**", "*.index")) ]
    selector = OutputSelector(outputs = fns)
    return selector
