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

class TransferFunction(object):
    def __init__(self, x_bounds, nbins=256):
        self.nbins = nbins
        self.x_bounds = x_bounds
        self.x = na.linspace(x_bounds[0], x_bounds[1], nbins).astype('float64')
        self.y = na.zeros(nbins, dtype='float64')

    def add_gaussian(self, location, width, height):
        vals = height * na.exp(-(self.x - location)**2.0/width)
        self.y = na.clip(na.maximum(vals, self.y), 0.0, 1.0)

    def add_line(self, start, stop):
        x0, y0 = start
        x1, y1 = stop
        slope = (y1-y0)/(x1-x0)
        vals = na.zeros(self.x.shape, 'float64')
        vals[(self.x >= x0) & (self.x <= x1)] = \
            slope * (self.x - x0) + y0
        self.y = na.clip(na.maximum(vals, self.y), 0.0, 1.0)

    def plot(self, filename):
        import matplotlib;matplotlib.use("Agg");import pylab
        pylab.clf()
        pylab.plot(self.x, self.y, 'xk-')
        pylab.xlim(*self.x_bounds)
        pylab.ylim(0.0, 1.0)
        pylab.savefig(filename)

class ColorTransferFunction(object):
    def __init__(self, x_bounds, nbins=256):
        self.x_bounds = x_bounds
        self.nbins = nbins
        self.red = TransferFunction(x_bounds, nbins)
        self.green = TransferFunction(x_bounds, nbins)
        self.blue = TransferFunction(x_bounds, nbins)
        self.alpha = TransferFunction(x_bounds, nbins)
        self.funcs = (self.red, self.green, self.blue, self.alpha)
        self.light_dir = (0.3,-0.2,0.5)
        self.light_color = (0.10, 0.10, 0.10)
        self.use_light = 1

    def add_gaussian(self, location, width, height):
        for tf, v in zip(self.funcs, height):
            tf.add_gaussian(location, width, v)

    def plot(self, filename):
        import matplotlib;matplotlib.use("Agg");import pylab
        pylab.clf()
        for c,tf in zip(['r','g','b'], self.funcs):
            pylab.plot(tf.x, tf.y, '-' + c)
            pylab.fill(tf.x, tf.y, c, alpha=0.2)
        pylab.plot(self.alpha.x, self.alpha.y, '-k')
        pylab.fill(self.alpha.x, self.alpha.y, 'k', alpha=0.1)
        pylab.xlim(*self.x_bounds)
        pylab.ylim(0.0, 1.0)
        pylab.xlabel("Value")
        pylab.ylabel("Transmission")
        pylab.savefig(filename)

if __name__ == "__main__":
    tf = ColorTransferFunction((-20, -5))
    tf.add_gaussian(-16.0, 0.4, [0.2, 0.3, 0.1])
    tf.add_gaussian(-14.0, 0.8, [0.4, 0.1, 0.2])
    tf.add_gaussian(-10.0, 1.0, [0.0, 0.0, 1.0])
    tf.plot("tf.png")
