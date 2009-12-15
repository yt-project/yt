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

from TransferFunction import ColorTransferFunction

from enthought.traits.api import \
        HasTraits, Float, List, Instance, Button, Array, CArray, Range, \
        DelegatesTo, Property, Any, Code
from enthought.traits.ui.api import \
    View, Item, HSplit, VSplit, ListEditor, InstanceEditor, ValueEditor, \
    HGroup, VGroup, CodeEditor, TextEditor
from enthought.chaco.api import Plot, ArrayPlotData
from enthought.enable.component_editor import ComponentEditor

class TFGaussian(HasTraits):
    center = Range(low = 'left_edge',
                   high = 'right_edge')
    left_edge = DelegatesTo('tf')
    right_edge = DelegatesTo('tf')

    tf = Any
    
    width = Float(0.0)

    red = Range(0.0, 1.0, 0.2)
    green = Range(0.0, 1.0, 0.2)
    blue = Range(0.0, 1.0, 0.2)
    alpha = Range(0.0, 1.0, 0.2)

    traits_view = View(VGroup(
                         HGroup(
                    Item('center'),
                    Item('width')
                               ),
                         HGroup(
                    Item('red'),
                    Item('green'),
                    Item('blue'),
                    Item('alpha')
                               ),
                             ),
                       )

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

class TFColors(HasTraits):
    gaussians = List(Instance(TFGaussian))
    transfer_function = Instance(ColorTransferFunction)

    left_edge = Float
    right_edge = Float

    add_gaussian = Button
    generate_code = Button

    plot_data = Instance(ArrayPlotData)
    plot = Instance(Plot)

    traits_view = View(VGroup(
                         Item('plot', editor=ComponentEditor(),
                                      show_label=False, resizable=True),
                         Item("gaussians", style='custom',
                              editor=ListEditor(style='custom'),
                              show_label=False,
                             ),
                         HGroup(Item("left_edge"), Item("right_edge")),
                         HGroup(Item("add_gaussian", show_label = False),
                                Item("generate_code", show_label = False),
                                ),
                        ),
                       width=500, height=500,
                       resizable=True)

    def _plot_data_default(self):
        return ArrayPlotData(rx = (0.0, 1.0), ry = (0.0, 0.0),
                             gx = (0.0, 1.0), gy = (0.0, 0.0),
                             bx = (0.0, 1.0), by = (0.0, 0.0),
                             ax = (0.0, 1.0), ay = (0.0, 0.0),
                             lx = (0.0, 1.0), ly = (0.0, 0.0),
                             ux = (0.0, 1.0), uy = (1.0, 1.0))

    def _plot_default(self):
        p = Plot(self.plot_data)
        p.plot( ("rx", "ry"), type='line', color='red')
        p.plot( ("gx", "gy"), type='line', color='green')
        p.plot( ("bx", "by"), type='line', color='blue')
        p.plot( ("ax", "ay"), type='line', color='black')
        p.plot( ("lx", "ly"), type='line', color='black')
        p.plot( ("ux", "uy"), type='line', color='black')
        return p

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

    def _generate_code_fired(self):
        val = ""
        val += "from yt.extensions.volume_rendering import ColorTransferFunction\n"
        val += "\n"
        val += "tf = ColorTransferFunction((%s, %s))\n" % (self.left_edge, self.right_edge)
        for g in self.gaussians:
            val += "tf.add_gaussian(%s, %s, (%s, %s, %s, %s))\n" % (
                    g.center, g.width, g.red, g.green, g.blue, g.alpha)
        val += "\n"
        snipped = CodeSnippet(snippet=val)
        snipped.edit_traits()

class CodeSnippet(HasTraits):
    snippet = Code

    # Not using code editor, because it doesn't work for me with Qt4
    traits_view = View(Item("snippet", editor=TextEditor(multi_line=True)),
                       width=800, height=600, resizable=True)

if __name__ == "__main__":
    tfc = TFColors()
    tfc.configure_traits()
