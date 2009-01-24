"""
The slowest ray-tracer known to man.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *

class SlowRayTracer(object):
    def __init__(self, pf, axis, bounds, res, func, reverse = False):
        self.pf = pf
        self.axis = axis
        self.x_bounds, self.y_bounds = bounds
        self.res = res
        self.func = func
        self.reverse = reverse
        self.buff = na.zeros(res, dtype='float64')
        self.fill_buff()

    def _trace_single_ray(self, x, y):
        ray = self.pf.h.ortho_ray(self.axis, (x, y))
        order = ray[axis_names[self.axis]].argsort()
        if self.reverse: order = order[::-1]
        return self.func(ray, order)

    def fill_buff(self):
        #pb = get_pbar("Ray tracing ", na.prod(self.res))
        xs = na.linspace(self.x_bounds[0], self.x_bounds[1], self.res[0])
        ys = na.linspace(self.y_bounds[0], self.y_bounds[1], self.res[1])
        #counter = 0
        for i, x in enumerate(xs.ravel()):
            #counter += 1
            print i
            for j, y in enumerate(ys.ravel()):
                #pb.update(counter)
                self.buff[i,j] = self._trace_single_ray(x,y)
        #pb.finish()

    def __getitem__(self, item):
        return self.buff[item]

if __name__ == "__main__":
    def raw_proj(ray, order):
        return na.sum((ray["Density"]*ray["dx"])[order])
        #return na.sum((ray["Density"]*ray["dx"])[order])
    import yt.logger
    yt.logger.disable_stream_logging()
    pf = EnzoStaticOutput("/Users/matthewturk/Research/data/galaxy1200.dir/galaxy1200")
    srt = SlowRayTracer(pf, 0, ((0.0125,1.0-0.0125),(0.0125,1.0-0.0125)),
                        (256, 256), raw_proj, reverse=False)
    import pylab
    pylab.imshow(na.log10(srt.buff), interpolation='nearest', origin='lower')
    pylab.savefig("srt_raw_x.png")
    pylab.clf()
