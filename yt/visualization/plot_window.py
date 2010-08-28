"""
A plotting mechanism based on the idea of a "window" into the data.

Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 J. S. Oishi.  All Rights Reserved.

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

import color_maps
from image_writer import \
    write_image
from fixed_resolution import \
    FixedResolutionBuffer
import matplotlib.pyplot

def invalidate_data(f):
    def newfunc(*args, **kwargs):
        args[0]._data_valid = False
        return f(*args, **kwargs)
    return newfunc

def invalidate_plot(f):
    def newfunc(*args, **kwargs):
        args[0]._plot_valid = False
        return f(*args, **kwargs)
    return newfunc

class PlotWindow(object):
    def __init__(self, data_source, bounds, buff_size=(800,800), antialias = True):
        r"""
        PlotWindow(data_source, bounds, buff_size=(800,800), antialias = True)
        
        A ploting mechanism based around the concept of a window into a
        data source. It can have arbitrary fields, each of which will be
        centered on the same viewpoint, but will have individual zlimits. 
        
        The data and plot are updated separately, and each can be
        invalidated as the object is modified.
        
        Data is handled by a FixedResolutionBuffer object.
        """
        self.data_source = data_source
        self.bounds = bounds
        self.buff_size = buff_size
        self.antialias = True

        self.plots = {}
        self._recreate_frb()
        self._frb._get_data_source_fields()
        self._setup_plots()
        self._data_valid = True

    def __getitem__(self, item):
        return self.plots[item]

    def _recreate_frb(self):
        try:
            self._frb = FixedResolutionBuffer(self.data_source, 
                                              self.bounds, self.buff_size, 
                                              self.antialias)
        except:
            raise RuntimeError("Failed to repixelize.")
        self._data_valid = True

    def _setup_plots(self):
        for f in self.fields:
            self.plots[f] = YtWindowPlot(self._frb[f])

    @property
    def fields(self):
        return self._frb.data.keys()

    def save(self):
        for p in self.plots:
            self.plots[p].save("blah")

    # def save(self):
    #     """
    #     REPLACE THIS
    #     """
    #     for field in self._frb.data.keys():
    #         name = "%s.png" % field
    #         print "writing %s" % name
    #         write_image(self._frb[field],name)
    
    @invalidate_data
    def pan(self):
        pass

    @invalidate_plot
    def set_cmap(self):
        pass

    @invalidate_data
    def set_field(self):
        pass

    @invalidate_data
    def set_window(self):
        pass

    @invalidate_data
    def set_width(self):
        pass

    # @invalidate_plot
    # def set_zlim(self):
    #     pass

    @invalidate_data
    def set_antialias(self,aa):
        self.antialias = aa

class YtPlot(object):
    """A base class for all yt plots. It should abstract the actual
    plotting engine completely, allowing plotting even without matplotlib. 

    YtPlot and the classes that derive from it are *by design* limited
    and designed for rapid, production quality plot production, rather
    than full generality. If you require more customization of the end
    result, these objects are designed to return to you the basic data
    so you the user can insert them into a matplotlib figure on your
    own outside of the YtPlot class.

    """
    axis_names = {}
    datalabel = None
    figure = None
    def __init__(self, field, size=(10,8)):
        self.__setup_from_field(field)
        self._plot_valid = True
        self.figure = matplotlib.pyplot.figure(figsize=size)
        self.axes = self.figure.add_subplot(1,1,1)

    def __setup_from_field(self, field):
        #self.set_log_field(self.pf.field_info[field].take_log)
        self.axis_names["X"] = None
        self.axis_names["Y"] = None
        self.axis_names["Z"] = field

    def save(self,name):
        print "saving plot %s" % name
        self.figure.savefig('%s.png' % name)

class Yt2DPlot(YtPlot):
    zmin = None
    zmax = None
    cmap = 'algae'
    zlabel = None

    # def __init__(self, data):
    #     pass

    @invalidate_plot
    def set_zlim(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax

class YtWindowPlot(Yt2DPlot):
    def __init__(self, data, size=(10,8)):
        YtPlot.__init__(self, data, size)
        self.__init_image(data)

    def __init_image(self, data):
        self.image = self.axes.imshow(data,cmap=self.cmap)

class YtProfilePlot(Yt2DPlot):
    def __init__(self):
        pass
