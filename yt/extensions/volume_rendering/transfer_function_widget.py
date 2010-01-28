"""
Simple transfer function editor

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

import numpy as na
import cPickle
from TransferFunction import ColorTransferFunction

from enthought.traits.api import \
        HasTraits, Float, List, Instance, Button, Array, CArray, Range, \
        DelegatesTo, Property, Any, Code, Callable
from enthought.traits.ui.api import \
    View, Item, HSplit, VSplit, ListEditor, InstanceEditor, ValueEditor, \
    HGroup, VGroup, CodeEditor, TextEditor
from enthought.chaco.api import Plot, ArrayPlotData
from enthought.enable.component_editor import ComponentEditor
import enthought.pyface.api as pyface

class TFGaussian(HasTraits):
    center = Range(low = 'left_edge',
                   high = 'right_edge')
    left_edge = DelegatesTo('tf')
    right_edge = DelegatesTo('tf')

    tf = Any
    
    width = Property
    rwidth = Range(0.0, 0.5, 0.05)

    red = Range(0.0, 1.0, 0.5)
    green = Range(0.0, 1.0, 0.5)
    blue = Range(0.0, 1.0, 0.5)
    alpha = Range(0.0, 1.0, 1.0)

    traits_view = View(VGroup(
                         HGroup(
                    Item('center'),
                    Item('rwidth', label='Width')
                               ),
                         HGroup(
                    Item('red'),
                    Item('green'),
                    Item('blue'),
                    Item('alpha')
                               ),
                             ),
                       )

    def _get_width(self):
        width = self.rwidth * (self.tf.right_edge - self.tf.left_edge)
        return width

    def _center_default(self):
        return (self.left_edge + self.right_edge)/2.0

    def _width_default(self):
        return (self.right_edge - self.left_edge)/20.0

    def _red_changed(self):
        self.tf._redraw()

    def _green_changed(self):
        self.tf._redraw()

    def _blue_changed(self):
        self.tf._redraw()

    def _alpha_changed(self):
        self.tf._redraw()

    def _center_changed(self):
        self.tf._redraw()

    def _height_changed(self):
        self.tf._redraw()

    def _rwidth_changed(self):
        self.tf._redraw()

class TFColors(HasTraits):
    gaussians = List(Instance(TFGaussian))
    transfer_function = Instance(ColorTransferFunction)

    left_edge = Float(0.0)
    right_edge = Float(10.0)

    add_gaussian = Button
    run_routine = Button
    save_function = Button

    routine = Callable

    plot_data = Instance(ArrayPlotData)
    image_data = Instance(ArrayPlotData)
    vr_image_data = Instance(ArrayPlotData)

    plot = Instance(Plot)
    image_plot = Instance(Plot)
    vr_image_plot = Instance(Plot)

    traits_view = View(VGroup(
                         HGroup(
                       VGroup(
                         Item('image_plot', editor=ComponentEditor(),
                                      show_label=False, resizable=True),
                         Item('plot', editor=ComponentEditor(),
                                      show_label=False, resizable=True),
                             ),
                         Item('vr_image_plot', editor=ComponentEditor(),
                                      show_label=False, resizable=True,
                                      width=512, height=512)),
                         Item("gaussians", style='custom',
                              editor=ListEditor(style='custom'),
                              show_label=False,
                             ),
                         HGroup(Item("left_edge"), Item("right_edge")),
                         HGroup(Item("add_gaussian", show_label = False),
                                Item("run_routine", show_label = False),
                                Item("save_function", show_label = False),
                                ),
                        ),
                       width=960, height=800,
                       resizable=True)

    def _plot_data_default(self):
        return ArrayPlotData(rx = (0.0, 1.0), ry = (0.0, 0.0),
                             gx = (0.0, 1.0), gy = (0.0, 0.0),
                             bx = (0.0, 1.0), by = (0.0, 0.0),
                             ax = (0.0, 1.0), ay = (0.0, 0.0),
                             lx = (0.0, 1.0), ly = (0.0, 0.0),
                             ux = (0.0, 1.0), uy = (1.0, 1.0))

    def _image_data_default(self):
        return ArrayPlotData(image_data = na.zeros((40,256,4), dtype='uint8'))

    def _vr_image_data_default(self):
        return ArrayPlotData(vr_image_data = na.zeros((512,512,3), dtype='uint8'))

    def _plot_default(self):
        p = Plot(self.plot_data)
        p.plot( ("rx", "ry"), type='line', color='red')
        p.plot( ("gx", "gy"), type='line', color='green')
        p.plot( ("bx", "by"), type='line', color='blue')
        p.plot( ("ax", "ay"), type='line', color='black')
        p.plot( ("lx", "ly"), type='line', color='black')
        p.plot( ("ux", "uy"), type='line', color='black')
        return p

    def _image_plot_default(self):
        plot = Plot(self.image_data, default_origin="top left")
        #plot.x_axis.orientation = "top"
        img_plot = plot.img_plot("image_data")[0]

        plot.bgcolor = "black"
        return plot

    def _vr_image_plot_default(self):
        plot = Plot(self.vr_image_data, default_origin="top left")
        #plot.x_axis.orientation = "top"
        img_plot = plot.img_plot("vr_image_data")[0]

        plot.bgcolor = "black"
        return plot

    def _add_gaussian_fired(self):
        self.gaussians.append(TFGaussian(tf = self))

    def _redraw(self):
        self.transfer_function = ColorTransferFunction(
                                  (self.left_edge, self.right_edge))
        for g in self.gaussians:
            self.transfer_function.add_gaussian(g.center, g.width,
                                             (g.red, g.green, g.blue, g.alpha))
        for f, c in zip(self.transfer_function.funcs, "rgba"):
            self.plot_data["%sx" % c] = f.x
            self.plot_data["%sy" % c] = f.y

        # Now we update the image describing the colors
        # This makes the assumption that all the x values are the same
        image = na.zeros((40, self.transfer_function.nbins, 4), dtype='uint8')
        for i,f in enumerate(self.transfer_function.funcs):
            image[:,:,i] = (f.y[None,:] * 255).astype('uint8')
        self.image_data["image_data"] = image

    def _run_routine_fired(self):
        img_data = self.routine(self.transfer_function)
        self.vr_image_data['vr_image_data'] = img_data

    def _save_function_fired(self):
        self._redraw()
        dlg = pyface.FileDialog(
            action='save as',
            wildcard="*.ctf",
        )
        if dlg.open() == pyface.OK:
            print "Saving:", dlg.path
            tf = self.transfer_function
            f = open(dlg.path, "wb")
            cPickle.dump(tf, f)

if __name__ == "__main__":
    tfc = TFColors()
    tfc.configure_traits()
