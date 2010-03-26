"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

# Major library imports
from vm_panner import VariableMeshPanner
from numpy import linspace, meshgrid, pi, sin, mgrid

# Enthought library imports
from enthought.enable.api import Component, ComponentEditor, Window
from enthought.traits.api import HasTraits, Instance, Button, Any
from enthought.traits.ui.api import Item, Group, View

# Chaco imports
from enthought.chaco.api import ArrayPlotData, jet, Plot
from enthought.chaco.image_data import FunctionImageData
from enthought.chaco.tools.api import PanTool, ZoomTool
from enthought.chaco.tools.image_inspector_tool import ImageInspectorTool, \
     ImageInspectorOverlay
from enthought.chaco.function_data_source import FunctionDataSource

class ImagePixelizerHelper(object):
    index = None
    def __init__(self, panner):
        self.panner = panner

    def __call__(self, low, high):
        b = self.panner.set_low_high(low, high)
        if self.index is not None:
            num_x_ticks = b.shape[0] + 1
            num_y_ticks = b.shape[1] + 1
            xs = mgrid[low[0]:high[0]:num_x_ticks*1j]
            ys = mgrid[low[1]:high[1]:num_y_ticks*1j]
            self.index.set_data( xs, ys )
        return b

class VariableMeshPannerView(HasTraits):

    plot = Instance(Plot)
    pd = Instance(ArrayPlotData)
    panner = Instance(VariableMeshPanner)
    fid = Instance(FunctionImageData)
    limits = Button
    helper = Any
    
    traits_view = View(
                    Group(
                        Item('plot', editor=ComponentEditor(size=(512,512)), 
                             show_label=False),
                        Item('limits', show_label=False),
                        orientation = "vertical"),
                    width = 800, height=800,
                    resizable=True, title="Pan and Scan",
                    )
    
    def __init__(self, **kwargs):
        super(VariableMeshPannerView, self).__init__(**kwargs)
        # Create the plot
        pd = ArrayPlotData()
        plot = Plot(pd)
        self.pd = pd
        helper = ImagePixelizerHelper(self.panner)
        fid = FunctionImageData(func = helper)
        fid._data = self.panner.buffer
        self.fid = fid
        bounds = self.panner.bounds
        pd.set_data("imagedata", fid)

        img_plot = plot.img_plot("imagedata", colormap=jet,
                                 interpolation='nearest',
                                 xbounds=(0.0, 1.0),
                                 ybounds=(0.0, 1.0))[0]
        helper.index = img_plot.index
        self.helper = helper

        fid.data_range = plot.range2d

        plot.tools.append(PanTool(img_plot))
        zoom = ZoomTool(component=img_plot, tool_mode="box", always_on=False)
        plot.overlays.append(zoom)
        imgtool = ImageInspectorTool(img_plot)
        img_plot.tools.append(imgtool)
        overlay = ImageInspectorOverlay(component=img_plot, image_inspector=imgtool,
                                        bgcolor="white", border_visible=True)

        img_plot.overlays.append(overlay)
        self.plot = plot

    def _limits_fired(self):
        print self.pd["imagedata"].min(), self.pd["imagedata"].max(),
        print self.fid.data.min(), self.fid.data.max()
