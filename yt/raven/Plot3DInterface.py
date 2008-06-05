"""
This is an interface to 
U{S2PLOT <http://http://astronomy.swin.edu.au/s2plot/index.php?title=S2PLOT>}
to plot uniform-spaced grids, derived from AMR data.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

from yt.raven import *
import s2plot

class S2PlotNotInitialized(Exception):
    pass

def must_have_s2plot(func):
    def check_started(obj, *args, **kwargs):
        if not obj.started:
            raise S2PlotNotInitialized()
        func(obj, *args, **kwargs)
    return check_started

class VolumeRenderingCube(object):
    def __init__(self, pf, center=None, width=1, unit='1',
                 field='Density', dims=128, take_log=True,
                 window_opts="/S2MONO", cmap="rainbow",
                 amin=0.01, amax=0.1):
        self.pf = pf
        self.width = width/pf[unit]
        if center is None: center = pf.h.find_max("Density")[1]
        self.center = center
        self.field = field
        self.dims = dims
        dx = self.width / dims
        self.max_level = na.unique(pf.h.gridDxs[pf.h.gridDxs>=dx]).argmax()
        self.__setup_data(take_log)
        self.vrid = None
        self.isoids = []
        self.window_opts = window_opts
        self.cmap = cmap
        self._amin, self._amax = amin, amax
        self.tr = na.array([ 0.0, (1.0-0.0)/(self.dims-1.0), 0, 0,
                             0.0, 0, (1.0-0.0)/(self.dims-1.0), 0,
                             0.0, 0, 0, (1.0-0.0)/(self.dims-1.0) ])
        self.started = False
        
    def __setup_data(self, take_log):
        self.data_grid = self.pf.h.covering_grid(self.max_level,
            self.center - self.width/2.0,
            self.center + self.width/2.0,
            [self.dims]*3, fields=[self.field])
        if take_log: self.data=na.log10(self.data_grid[self.field])
        else: self.data=self.data_grid[self.field]
        self._dmin = self.data.min()
        self._dmax = self.data.max()

    @must_have_s2plot
    def add_isosurfaces(self, number=None, vals=None,
                        log_space=True, cmap="jet",
                        amin=None, amax=None):
        if amin is None: amin = self._amin
        if amax is None: amax = self._amax
        cm = matplotlib.cm.get_cmap(cmap)
        if number is None and val is None:
            raise ValueError("You have to supply either number or vals")
        if number is None: number = len(vals)
        if vals is None:
            if log_space: func=na.logspace
            else: func=na.linspace
            vals = func(self._dmin, self._dmax, number)
            if log_space: vals = na.log10(vals)
        for val,a in zip(vals, na.linspace(amin, amax, number)):
            self.isoids.append(
                self.__add_isosurface(val, a, cm))
        for id in self.isoids: s2plot.ns2dis(id,0)

    def __add_isosurface(self, val, alpha, cm):
        scaled_val = ((val-self._dmin)/(self._dmax-self._dmin))
        r,g,b,a = cm(scaled_val)
        nx = self.dims-1
        print "Adding",val,scaled_val,alpha,r,g,b
        return s2plot.ns2cis(self.data, self.dims, self.dims, self.dims,
                             0, nx, 0, nx, 0, nx, self.tr, val,
                             1, 't', alpha, r,g,b)
                            
    def run(self):
        self.__setup_s2plot()
        self.__setup_volrendering()
        self.__register_callbacks()
        self.__start_rendering()

    def restart(self):
        self.__start_rendering()

    def __setup_s2plot(self):
        dx = 1.0/(self.dims-1.0)
        s2plot.s2opendo(self.window_opts)
        s2plot.s2swin(-dx,1.0+dx,-dx,1.0+dx,-dx,1.0+dx) # We mandate cube
        s2plot.s2box("BCDE",0,0,"BCDE",0,0,"BCDE",0,0)
        s2plot.s2scir(1000,2000)            # Set colour range
        s2plot.s2icm(self.cmap,1000,2000)   # Install colour map
        amb = {'r':0.8, 'g':0.8, 'b':0.8}   # ambient light
        s2plot.ss2srm(s2plot.SHADE_FLAT);   # Set shading type to FLAT
        s2plot.ss2sl(amb, 0, None, None, 0) # Ambient lighting only
        self.started = True

    def __setup_volrendering(self):
        nd = self.dims
        self.vrid = s2plot.ns2cvr(self.data, nd, nd, nd,
                           0, nd-1, 0, nd-1, 0, nd-1, 
                           self.tr, 't',
                           self._dmin, self._dmax, self._amin, self._amax)
        
    def __register_callbacks(self):
        # More should go here for functional changes to the object
        s2plot.cs2scb(self.__my_callback)                           # Install a dynamic callback

    def __start_rendering(self):
        s2plot.s2disp(-1, 1)

    def __my_callback(self, t, kc):
        s2plot.ds2dvr(self.vrid, 0)
