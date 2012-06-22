"""
These widget objects provide the interaction to yt tasks.

Author: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
import weakref

class RenderingScene(object):
    _camera = None
    _tf = None
    
    def __init__(self, pf, camera=None, tf=None):
        self.pf = weakref.proxy(pf)
        self._camera = camera
        self._tf = tf

        self.center = self.pf.domain_center
        self.normal_vector = na.array([0.7,1.0,0.3])
        self.north_vector = [0.,0.,1.]
        self.steady_north = True
        self.fields = ['Density']
        self.log_fields = [True]
        self.l_max = 0
        self.mi = None
        self.ma = None
        if self._tf is None:
            self._new_tf(self.pf)

        if self._camera is None:
            self._new_camera(self.pf)

    def _new_tf(self, pf, mi=None, ma=None, nbins=1024):
        if mi is None or ma is None:
            roi = self.pf.h.region(self.center, self.center-self.width, self.center+self.width)
            self.mi, self.ma = roi.quantities['Extrema'](self.fields[0])[0]
            if self.log_fields[0]:
                self.mi, self.ma = na.log10(self.mi), na.log10(self.ma)

        self._tf = ColorTransferFunction((self.mi-2, self.ma+2), nbins=nbins)

    def add_contours(self, n_contours=7, contour_width=0.05, colormap='kamae'):
        self._tf.add_layers(n_contours=n_contours ,w=contour_width,
                                  col_bounds = (self.mi,self.ma), 
                                  colormap=colormap)

    def _new_camera(self, pf):
        del self._camera
        self._camera = self.pf.camera(self.center, self.normal_vector, 
                                      self.width, self.resolution, self._tf,
                                      north_vector=self.north_vector,
                                      steady_north=self.steady_north,
                                      fields=self.fields, log_fields=self.log_fields,
                                      l_max=self.l_max)
    def snapshot(self):
        return self._camera.snapshot()


def get_corners(pf, max_level=Non):
    corners = pf.h.grid_corners.tolist()
    levels = pf.h.grid_levels.to_list()
    return corners, levels    

def get_grid_vertices(pf, max_level=None)
    if max_level is None: max_level = pf.h.max_level
    fns = []
    for lvl in range(max_level):
        fn = "%s_lvl_%02i_grids.vtk"%(pf, lvl)
        f = open(fn, "w")

        f.write("""# vtk DataFile Version 3.0
        vtk output
        ASCII
        DATASET POLYDATA
        POINTS %s float
        """ % (pf.h.num_grids * 8))

        for gi in xrange(pf.h.num_grids):
            if pf.h.grid_levels[gi][0] != lvl: continue
            gc = pf.h.grid_corners[:, :, gi] * 100.0
            for i in range(8):
                f.write("%0.9f %0.9f %0.9f\n" % (gc[i, 0], gc[i, 1], gc[i, 2]))

        f.write("LINES %s %s\n" % (pf.h.num_grids * 12, 3 * pf.h.num_grids * 12))

        order1 = (0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3)
        order2 = (1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7)

        for gi in xrange(pf.h.num_grids):
            if pf.h.grid_levels[gi][0] != lvl: continue
            offset = 8 * gi
            for p1, p2 in zip(order1, order2):
                f.write("2 %s %s\n" % (p1 + offset, p2 + offset))
        f.close()
        fns.append(fn)
    return fn

def get_isocontour(pf, field, value=None):

    dd = pf.h.all_data()
    if value is None:
        mi, ma = na.log10(dd.quantities["Extrema"]("Density")[0])
        value = 10.**((mi + ma)/2.)
    vert = dd.extract_isocontours("Density", value)
    na.multiply(vert, 100, vert)

    fn = "%s_surface.vtk" % pf
    f = open(fn, "w")
    f.write("vtk output\nASCII\nDATASET POLYDATA\nPOINTS %s float\n" % (vert.shape[0]))
    for v in vert:
        f.write("%0.14f %0.14f %0.14f\n" % (v[0], v[1], v[2]))
    f.write("VERTICES %s %s\n" % (vert.shape[0]/3, vert.shape[0] + vert.shape[0]/3))
    for i in range(vert.shape[0]/3):
        f.write("3 %i %i %i\n" % (i*3, i*3+1, i*3+2))
    f.close()
    return fn

def get_streamlines(pf):
    from yt.visualization.api import Streamlines
    streamlines = Streamlines(pf, pf.domain_center) 
    streamlines.integrate_through_volume()
    stream = streamlines.path(0)
    matplotlib.pylab.semilogy(stream['t'], stream['Density'], '-x')


