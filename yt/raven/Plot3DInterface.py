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

class VolumeRendering(object):
    xyz = None
    def __init__(self, data, take_log=True,
                 window_opts="/S2MONO", cmap="rainbow",
                 amin=0.0, amax=0.1, bounds = None):
        if bounds is None:
            self.x0,self.x1, self.y0,self.y1, self.z0,self.z1 = \
                0.0,1.0, 0.0,1.0, 0.0,1.0
        else:
            self.x0,self.x1, self.y0,self.y1, self.z0,self.z1 = \
                bounds
        self.__setup_data(data, take_log)
        self.vrid = None
        self.isoids = []
        self.window_opts = window_opts
        self.cmap = cmap
        self._amin, self._amax = amin, amax
        nx, ny, nz = self.dims
        self.tr = na.array([
            self.x0, (self.x1-self.x0)/(nx-1.0), 0               , 0,
            self.y0, 0               , (self.y1-self.y0)/(ny-1.0), 0,
            self.z0, 0               , 0               , (self.z1-self.z0)/(nz-1.0)
        ])
        self.started = False
        
    def __setup_data(self, data, take_log):
        self.dims = na.array(data.shape)
        if take_log: self.data=na.log10(data).copy()
        else: self.data=data.copy()
        self._dmin = self.data.min()
        if take_log: self._dmin=na.log10(data[data>0]).min()
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
            vals = func(self._dmin, self._dmax, number+2)[1:-1]
            if log_space: vals = na.log10(vals)
        for val,a in zip(vals, na.linspace(amin, amax, number)):
            self.isoids.append(
                self.__add_isosurface(val, a, cm))
        for id in self.isoids: s2plot.ns2dis(id,0)

    def __add_isosurface(self, val, alpha, cm):
        scaled_val = ((val-self._dmin)/(self._dmax-self._dmin))
        r,g,b,a = cm(scaled_val)
        nx,ny,nz = self.dims
        mylog.info("Adding isosurface at %0.5e (%0.3e) with alpha %0.2f (%0.2f %0.2f %0.2f)",
                   val,scaled_val,alpha, r,g,b)
        return s2plot.ns2cis(self.data, nx, ny, nz,
                             0, nx-1, 0, ny-1, 0, ny-1,
                             self.tr, val,
                             1, 's', alpha, r,g,b)
                            
    def run(self, pre_call=None):
        self.__setup_s2plot()
        self.__setup_volrendering()
        self.__register_callbacks()
        self._setup_labels()
        if pre_call is not None: pre_call(self)
        s2plot.s2disp(-1, 1)

    def restart(self):
        s2plot.s2disp(-1, 0)

    def _setup_labels(self):
        pass

    def __setup_s2plot(self):
        dx,dy,dz = 1.0/(self.dims-1.0)
        s2plot.s2opendo(self.window_opts)
        s2plot.s2swin(-dx/1.0,1.0+dx/1.0,
                      -dy/1.0,1.0+dy/1.0,
                      -dz/1.0,1.0+dz/1.0)
        #opts = "BCDETMNOPQ"
        opts = "BCDE"
        s2plot.s2box(opts,0,0,opts,0,0,opts,0,0)
        s2plot.s2scir(1000,2000)            # Set colour range
        s2plot.s2icm(self.cmap,1000,2000)   # Install colour map
        amb = {'r':0.8, 'g':0.8, 'b':0.8}   # ambient light
        s2plot.ss2srm(s2plot.SHADE_FLAT);   # Set shading type to FLAT
        s2plot.ss2sl(amb, 0, None, None, 0) # Ambient lighting only
        self.started = True

    def __setup_volrendering(self):
        ndx,ndy,ndz = self.dims
        self.vrid = s2plot.ns2cvr(self.data, ndx, ndy, ndz,
                           0, ndx-1, 0, ndy-1, 0, ndz-1, 
                           self.tr, 's',
                           self._dmin, self._dmax, self._amin, self._amax)
        
    def __register_callbacks(self):
        # More should go here for functional changes to the object
        s2plot.cs2scb(self.__my_callback)                           # Install a dynamic callback

    def __my_callback(self, t, kc):
        s2plot.ds2dvr(self.vrid, 0)
        if self.xyz is not None and kc % 2 == 0 and kc > 0:
            s2plot.s2sci(s2plot.S2_PG_LTGREY)
            s2plot.s2pt(*self.xyz)

class VolumeRenderingDataCube(VolumeRendering):
    def __init__(self, pf, center=None, width=1, unit='1',
                 field='Density', dims=128, smooth_data=True,
                 **kwargs):
        self.pf = pf
        self.width = width/pf[unit]
        if center is None: center = pf.h.find_max("Density")[1]
        self.center = center
        self._use_smoothed = smooth_data
        self.field = field
        self.dims = dims
        dx = self.width / dims
        self.max_level = na.unique(pf.h.gridDxs[pf.h.gridDxs>=dx]).argmax()+1
        self.data_grid = self.__get_data()
        self.xyz = None
        VolumeRendering.__init__(self, self.data_grid[field], **kwargs)
        
    def __get_data(self):
        if self._use_smoothed: cl = self.pf.h.smoothed_covering_grid
        else: cl = self.pf.h.covering_grid
        data_grid = cl(
            level=self.max_level,
            left_edge=self.center - self.width/2.0,
            right_edge=self.center + self.width/2.0,
            dims=[self.dims]*3, fields=[self.field])
        return data_grid

class VolumeRendering3DProfile(VolumeRendering):
    def __init__(self, profile, field, **kwargs):
        self.profile = profile
        self.field = field
        if 'bounds' not in kwargs:
            kwargs['bounds'] = self._setup_bounds()
        VolumeRendering.__init__(self, self.profile[field], **kwargs)

    def _setup_bounds(self):
        nolog = lambda a: a
        self.bf = [[nolog, self.profile.x_bin_field], \
                   [nolog, self.profile.y_bin_field], \
                   [nolog, self.profile.z_bin_field]]
        if self.profile._x_log: self.bf[0][0] = na.log10
        if self.profile._y_log: self.bf[1][0] = na.log10
        if self.profile._z_log: self.bf[2][0] = na.log10
        xmi, ymi, zmi = [f(self.profile[fn]).min() for f,fn in self.bf]
        xma, yma, zma = [f(self.profile[fn]).max() for f,fn in self.bf]
        return (xmi,xma,ymi,yma,zmi,zma)

    def _setup_labels(self):
        s2plot.s2env(self.x0,self.x1, self.y0,self.y1, self.z0,self.z1, 0, 1)
        s2plot.s2lab(self.bf[0][1],self.bf[1][1],self.bf[2][1],self.field)

    def setup_plot_points(self):
        xyz = [f[0](self.profile._data_source[f[1]]) for f in self.bf]
        self.xyz = [xyz[0].size] + xyz +  [1]
