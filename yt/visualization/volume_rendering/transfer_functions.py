"""
Simple transfer function editor

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
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
from matplotlib.cm import get_cmap

from yt.funcs import *

from yt.utilities.physical_constants import *

class TransferFunction(object):
    def __init__(self, x_bounds, nbins=256):
        r"""A transfer function governs the transmission of emission and
        absorption through a volume.

        Transfer functions are defined by boundaries, bins, and the value that
        governs transmission through that bin.  This is scaled between 0 and 1.
        When integrating through a volume. the value through a given cell is
        defined by the value calculated in the transfer function.

        Parameters
        ----------
        x_bounds : tuple of floats
            The min and max for the transfer function.  Values below or above
            these values are discarded.
        nbins : int
            How many bins to calculate; in betwee, linear interpolation is
            used, so low values are typically fine.

        Notes
        -----
        Typically, raw transfer functions are not generated unless particular
        and specific control over the integration is desired.  Usually either
        color transfer functions, where the color values are calculated from
        color tables, or multivariate transfer functions are used.
        """
        self.pass_through = 0
        self.nbins = nbins
        self.x_bounds = x_bounds
        self.x = na.linspace(x_bounds[0], x_bounds[1], nbins).astype('float64')
        self.y = na.zeros(nbins, dtype='float64')
        self.grad_field = -1
        self.light_source_v = self.light_source_c = na.zeros(3, 'float64')

    def add_gaussian(self, location, width, height):
        r"""Add a Gaussian distribution to the transfer function.

        Typically, when rendering isocontours, a Guassian distribution is the
        easiest way to draw out features.  The spread provides a softness.
        The values are calculated as :math:`f(x) = h \exp{-(x-x_0)^2 / w}`.

        Parameters
        ----------
        location : float
            The centroid of the Gaussian (:math:`x_0` in the above equation.)
        width : float
            The relative width (:math:`w` in the above equation.)
        height : float
            The peak height (:math:`h` in the above equation.)  Note that while
            values greater 1.0 will be accepted, the values of the transmission
            function are clipped at 1.0.

        Examples
        --------

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, 1.0)
        """
        vals = height * na.exp(-(self.x - location)**2.0/width)
        self.y = na.clip(na.maximum(vals, self.y), 0.0, na.inf)

    def add_line(self, start, stop):
        r"""Add a line between two points to the transmission function.

        This will accept a starting point in (x,y) and an ending point in (x,y)
        and set the values of the transmission function between those x-values
        to be along the line connecting the y values.

        Parameters
        ----------
        start : tuple of floats
            (x0, y0), the starting point.  x0 is between the bounds of the
            transfer function and y0 must be between 0.0 and 1.0.
        stop : tuple of floats
            (x1, y1), the ending point.  x1 is between the bounds of the
            transfer function and y1 must be between 0.0 and 1.0.

        Examples
        --------
        This will set the transfer function to be linear from 0.0 to 1.0,
        across the bounds of the function.

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_line( (-10.0, 0.0), (-5.0, 1.0) )
        """
        x0, y0 = start
        x1, y1 = stop
        slope = (y1-y0)/(x1-x0)
        # We create a whole new set of values and then backout the ones that do
        # not satisfy our bounding box arguments
        vals = slope * (self.x - x0) + y0
        vals[~((self.x >= x0) & (self.x <= x1))] = 0.0
        self.y = na.clip(na.maximum(vals, self.y), 0.0, na.inf)

    def add_step(self, start, stop, value):
        r"""Adds a step function to the transfer function.

        This accepts a `start` and a `stop`, and then in between those points the
        transfer function is set to the maximum of the transfer function and
        the `value`.

        Parameters
        ----------
        start : float
            This is the beginning of the step function; must be within domain
            of the transfer function.
        stop : float
            This is the ending of the step function; must be within domain
            of the transfer function.
        value : float
            The value the transfer function will be set to between `start` and
            `stop`.  Note that the transfer function will *actually* be set to
            max(y, value) where y is the existing value of the transfer
            function.

        Examples
        --------
        Note that in this example, we have added a step function, but the
        Gaussian that already exists will "win" where it exceeds 0.5.

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-7.0, 0.01, 1.0)
        >>> tf.add_step(-8.0, -6.0, 0.5)
        """
        vals = na.zeros(self.x.shape, 'float64')
        vals[(self.x >= start) & (self.x <= stop)] = value
        self.y = na.clip(na.maximum(vals, self.y), 0.0, na.inf)

    def add_filtered_planck(self, wavelength, trans):
        vals = na.zeros(self.x.shape, 'float64')
        nu = clight/(wavelength*1e-8)
        nu = nu[::-1]

        for i,logT in enumerate(self.x):
            T = 10**logT
            # Black body at this nu, T
            Bnu = ((2.0 * hcgs * nu**3) / clight**2.0) / \
                    (na.exp(hcgs * nu / (kboltz * T)) - 1.0)
            # transmission
            f = Bnu * trans[::-1]
            # integrate transmission over nu
            vals[i] = na.trapz(f,nu)

        # normalize by total transmission over filter
        self.y = vals/trans.sum() #/na.trapz(trans[::-1],nu)
        #self.y = na.clip(na.maximum(vals, self.y), 0.0, 1.0)

    def plot(self, filename):
        r"""Save an image file of the transfer function.

        This function loads up matplotlib, plots the transfer function and saves.

        Parameters
        ----------
        filename : string
            The file to save out the plot as.

        Examples
        --------

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, 1.0)
        >>> tf.plot("sample.png")
        """
        import matplotlib;matplotlib.use("Agg");import pylab
        pylab.clf()
        pylab.plot(self.x, self.y, 'xk-')
        pylab.xlim(*self.x_bounds)
        pylab.ylim(0.0, 1.0)
        pylab.savefig(filename)

    def show(self):
        r"""Display an image of the transfer function

        This function loads up matplotlib and displays the current transfer function.

        Parameters
        ----------

        Examples
        --------

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, 1.0)
        >>> tf.show()
        """
        import matplotlib;import pylab
        pylab.clf()
        pylab.plot(self.x, self.y, 'xk-')
        pylab.xlim(*self.x_bounds)
        pylab.ylim(0.0, 1.0)
        pylab.draw()
        
    def clear(self):
        self.y[:]=0.0

class MultiVariateTransferFunction(object):
    def __init__(self):
        r"""This object constructs a set of field tables that allow for
        multiple field variables to control the integration through a volme.

        The integration through a volume typically only utilizes a single field
        variable (for instance, Density) to set up and control the values
        returned at the end of the integration.  For things like isocontours,
        this is fine.  However, more complicated schema are possible by using
        this object.  For instance, density-weighted emission that produces
        colors based on the temperature of the fluid.
        """
        self.n_field_tables = 0
        self.tables = [] # Tables are interpolation tables
        self.field_ids = [0] * 6 # This correlates fields with tables
        self.weight_field_ids = [-1] * 6 # This correlates 
        self.field_table_ids = [0] * 6
        self.weight_table_ids = [-1] * 6
        self.grad_field = -1
        self.light_source_v = self.light_source_c = na.zeros(3, 'float64')

    def add_field_table(self, table, field_id, weight_field_id = -1,
                        weight_table_id = -1):
        r"""This accepts a table describing integration.

        A "field table" is a tabulated set of values that govern the
        integration through a given field.  These are defined not only by the
        transmission coefficient, interpolated from the table itself, but the
        `field_id` that describes which of several fields the integration
        coefficient is to be calculated from.

        Parameters
        ----------
        table : `TransferFunction`
            The integration table to be added to the set of tables used during
            the integration.
        field_id : int
            Each volume has an associated set of fields.  This identifies which
            of those fields will be used to calculate the integration
            coefficient from this table.
        weight_field_id : int, optional
            If specified, the value of the field this identifies will be
            multiplied against the integration coefficient.
        weight_table_id : int, optional
            If specified, the value from the *table* this identifies will be
            multiplied against the integration coefficient.

        Notes
        -----
        This can be rather complicated.  It's recommended that if you are
        interested in manipulating this in detail that you examine the source
        code, specifically the function FIT_get_value in
        yt/_amr_utils/VolumeIntegrator.pyx.

        Examples
        --------
        This example shows how to link a new transfer function against field 0.
        Note that this by itself does not link a *channel* for integration
        against a field.  This is because the weighting system does not mandate
        that all tables contribute to a channel, only that they contribute a
        value which may be used by other field tables.

        >>> mv = MultiVariateTransferFunction()
        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian( -7.0, 0.01, 1.0)
        >>> mv.add_field_table(tf, 0)
        """
        self.tables.append(table)
        self.field_ids[self.n_field_tables] = field_id
        self.weight_field_ids[self.n_field_tables] = weight_field_id
        self.weight_table_ids[self.n_field_tables] = weight_table_id
        self.n_field_tables += 1

    def link_channels(self, table_id, channels = 0):
        r"""Link an image channel to a field table.

        Once a field table has been added, it can be linked against a channel (any
        one of the six -- red, green, blue, red absorption, green absorption, blue
        absorption) and then the value calculated for that field table will be
        added to the integration for that channel.  Not all tables must be linked
        against channels.

        Parameters
        ----------
        table_id : int
            The 0-indexed table to link.
        channels : int or list of ints
            The channel or channels to link with this table's calculated value.


        Examples
        --------
        This example shows how to link a new transfer function against field 0, and
        then link that table against all three RGB channels.  Typically an
        absorption (or 'alpha') channel is also linked.

        >>> mv = MultiVariateTransferFunction()
        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian( -7.0, 0.01, 1.0)
        >>> mv.add_field_table(tf, 0)
        >>> mv.link_channels(0, [0,1,2])
        """
        channels = ensure_list(channels)
        for c in channels:
            self.field_table_ids[c] = table_id

class ColorTransferFunction(MultiVariateTransferFunction):
    def __init__(self, x_bounds, nbins=256):
        r"""A complete set of transfer functions for standard color-mapping.

        This is the best and easiest way to set up volume rendering.  It
        creates field tables for all three colors, their alphas, and has
        support for sampling color maps and adding independent color values at
        all locations.  It will correctly set up the
        `MultiVariateTransferFunction`.

        Parameters
        ----------
        x_bounds : tuple of floats
            The min and max for the transfer function.  Values below or above
            these values are discarded.
        nbins : int
            How many bins to calculate; in betwee, linear interpolation is
            used, so low values are typically fine.
        """
        MultiVariateTransferFunction.__init__(self)
        self.x_bounds = x_bounds
        self.nbins = nbins
        # This is all compatibility and convenience.
        self.red = TransferFunction(x_bounds, nbins)
        self.green = TransferFunction(x_bounds, nbins)
        self.blue = TransferFunction(x_bounds, nbins)
        self.alpha = TransferFunction(x_bounds, nbins)
        self.funcs = (self.red, self.green, self.blue, self.alpha)

        # Now we do the multivariate stuff
        # We assign to Density, but do not weight
        for i,tf in enumerate(self.funcs[:3]):
            self.add_field_table(tf, 0, weight_table_id = 3)
            self.link_channels(i, i)
        self.add_field_table(self.funcs[3], 0)
        # We don't have a fifth table, so the value will *always* be zero.
        self.link_channels(4, [3,4,5])

    def add_gaussian(self, location, width, height):
        r"""Add a Gaussian distribution to the transfer function.

        Typically, when rendering isocontours, a Guassian distribution is the
        easiest way to draw out features.  The spread provides a softness.
        The values are calculated as :math:`f(x) = h \exp{-(x-x_0)^2 / w}`.

        Parameters
        ----------
        location : float
            The centroid of the Gaussian (:math:`x_0` in the above equation.)
        width : float
            The relative width (:math:`w` in the above equation.)
        height : list of 4 float
            The peak height (:math:`h` in the above equation.)  Note that while
            values greater 1.0 will be accepted, the values of the transmission
            function are clipped at 1.0.  This must be a list, and it is in the
            order of (red, green, blue, alpha).

        Examples
        --------
        This adds a red spike.

        >>> tf = ColorTransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, [1.0, 0.0, 0.0, 1.0])
        """
        alpha = height[3]
        for tf, v in zip(self.funcs, height):
            tf.add_gaussian(location, width, v)

    def add_step(self, start, stop, value):
        r"""Adds a step function to the transfer function.

        This accepts a `start` and a `stop`, and then in between those points the
        transfer function is set to the maximum of the transfer function and
        the `value`.

        Parameters
        ----------
        start : float
            This is the beginning of the step function; must be within domain
            of the transfer function.
        stop : float
            This is the ending of the step function; must be within domain
            of the transfer function.
        value : list of 4 floats
            The value the transfer function will be set to between `start` and
            `stop`.  Note that the transfer function will *actually* be set to
            max(y, value) where y is the existing value of the transfer
            function.  This must be a list, and it is in the order of (red,
            green, blue, alpha).


        Examples
        --------
        This adds a step function that will produce a white value at > -6.0.

        >>> tf = ColorTransferFunction( (-10.0, -5.0) )
        >>> tf.add_step(-6.0, -5.0, [1.0, 1.0, 1.0, 1.0])
        """
        for tf, v in zip(self.funcs, value):
            tf.add_step(start, stop, v)

    def plot(self, filename):
        r"""Save an image file of the transfer function.

        This function loads up matplotlib, plots all of the constituent
        transfer functions and saves.

        Parameters
        ----------
        filename : string
            The file to save out the plot as.

        Examples
        --------

        >>> tf = ColorTransferFunction( (-10.0, -5.0) )
        >>> tf.add_layers(8)
        >>> tf.plot("sample.png")
        """
        from matplotlib import pyplot
        from matplotlib.ticker import FuncFormatter
        pyplot.clf()
        ax = pyplot.axes()
        i_data = na.zeros((self.alpha.x.size, self.funcs[0].y.size, 3))
        i_data[:,:,0] = na.outer(na.ones(self.alpha.x.size), self.funcs[0].y)
        i_data[:,:,1] = na.outer(na.ones(self.alpha.x.size), self.funcs[1].y)
        i_data[:,:,2] = na.outer(na.ones(self.alpha.x.size), self.funcs[2].y)
        ax.imshow(i_data, origin='lower')
        ax.fill_between(na.arange(self.alpha.y.size), self.alpha.x.size * self.alpha.y, y2=self.alpha.x.size, color='white')
        ax.set_xlim(0, self.alpha.x.size)
        xticks = na.arange(na.ceil(self.alpha.x[0]), na.floor(self.alpha.x[-1]) + 1, 1) - self.alpha.x[0]
        xticks *= self.alpha.x.size / (self.alpha.x[-1] - self.alpha.x[0])
        ax.xaxis.set_ticks(xticks)
        def x_format(x, pos):
            return "%.1f" % (x * (self.alpha.x[-1] - self.alpha.x[0]) / (self.alpha.x.size) + self.alpha.x[0])
        ax.xaxis.set_major_formatter(FuncFormatter(x_format))
        yticks = na.linspace(0,1,5) * self.alpha.y.size
        ax.yaxis.set_ticks(yticks)
        def y_format(y, pos):
            return (y / self.alpha.y.size)
        ax.yaxis.set_major_formatter(FuncFormatter(y_format))
        ax.set_ylabel("Transmission")
        ax.set_xlabel("Value")
        pyplot.savefig(filename)

    def show(self, ax=None):
        r"""Display an image of the transfer function

        This function loads up matplotlib and displays the current transfer function.

        Parameters
        ----------

        Examples
        --------

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, 1.0)
        >>> tf.show()
        """
        from matplotlib import pyplot
        from matplotlib.ticker import FuncFormatter
        #pyplot.clf()
        if ax is None:
            ax = pyplot.axes()
        alpha = self.red.y + self.green.y + self.blue.y
        max_alpha = alpha.max()
        norm = max_alpha
        i_data = na.zeros((self.alpha.x.size, self.funcs[0].y.size, 3))
        i_data[:,:,0] = na.outer(na.ones(self.alpha.x.size), self.funcs[0].y/norm)
        i_data[:,:,1] = na.outer(na.ones(self.alpha.x.size), self.funcs[1].y/norm)
        i_data[:,:,2] = na.outer(na.ones(self.alpha.x.size), self.funcs[2].y/norm)
        ax.imshow(i_data, origin='lower', aspect='auto')

        #ax.fill_between(na.arange(self.alpha.y.size), alpha, y2=max_alpha, color='white')
        ax.plot(na.arange(self.alpha.y.size), alpha, 'w')
        ax.set_xlim(0, self.alpha.x.size)
        xticks = na.arange(na.ceil(self.alpha.x[0]), na.floor(self.alpha.x[-1]) + 1, 1) - self.alpha.x[0]
        xticks *= self.alpha.x.size / (self.alpha.x[-1] - self.alpha.x[0])
        if len(xticks) > 5:
            xticks = xticks[::len(xticks)/5]
        ax.xaxis.set_ticks(xticks)
        def x_format(x, pos):
            return "%.1f" % (x * (self.alpha.x[-1] - self.alpha.x[0]) / (self.alpha.x.size) + self.alpha.x[0])
        ax.xaxis.set_major_formatter(FuncFormatter(x_format))
        yticks = na.linspace(0,1,5) * max_alpha
        ax.yaxis.set_ticks(yticks)
        def y_format(y, pos):
            s = '%0.2f' % ( y )
            return s
        ax.yaxis.set_major_formatter(FuncFormatter(y_format))
        ax.set_ylim(0., max_alpha)
        ax.set_ylabel("Opacity")
        ax.set_xlabel("Value")

    def vert_cbar(self, ax=None, label=None):
        r"""Display an image of the transfer function

        This function loads up matplotlib and displays the current transfer function.

        Parameters
        ----------

        Examples
        --------

        >>> tf = TransferFunction( (-10.0, -5.0) )
        >>> tf.add_gaussian(-9.0, 0.01, 1.0)
        >>> tf.show()
        """
        from matplotlib import pyplot
        from matplotlib.ticker import FuncFormatter
        #pyplot.clf()
        if ax is None:
            ax = pyplot.axes()
        if label is None:
            label = ''
        alpha = self.red.y + self.green.y + self.blue.y
        max_alpha = alpha.max()
        norm = max_alpha
        i_data = na.zeros((self.alpha.x.size, self.funcs[0].y.size, 3))
        i_data[:,:,0] = na.outer(self.funcs[0].y/norm, na.ones(self.alpha.x.size))
        i_data[:,:,1] = na.outer(self.funcs[1].y/norm, na.ones(self.alpha.x.size))
        i_data[:,:,2] = na.outer(self.funcs[2].y/norm, na.ones(self.alpha.x.size))
        ax.imshow(i_data, origin='lower', aspect='auto')

        #ax.fill_between(na.arange(self.alpha.y.size), alpha, y2=max_alpha, color='white')
        #ax.plot(alpha, na.arange(self.alpha.y.size), 'w')
        ax.set_ylim(0, self.alpha.x.size)
        xticks = na.arange(na.ceil(self.alpha.x[0]), na.floor(self.alpha.x[-1]) + 1, 1) - self.alpha.x[0]
        xticks *= self.alpha.x.size / (self.alpha.x[-1] - self.alpha.x[0])
        if len(xticks) > 5:
            xticks = xticks[::len(xticks)/5]
        ax.yaxis.set_ticks(xticks)
        def x_format(x, pos):
            return "%.1f" % (x * (self.alpha.x[-1] - self.alpha.x[0]) / (self.alpha.x.size) + self.alpha.x[0])
        ax.yaxis.set_major_formatter(FuncFormatter(x_format))
        yticks = na.linspace(0,1,2,endpoint=True) * max_alpha
        ax.xaxis.set_ticks(yticks)
        #balls
        def y_format(y, pos):
            s = '%0.2f' % ( y )
            return s
        ax.xaxis.set_major_formatter(FuncFormatter(y_format))
        ax.set_xlim(0., max_alpha)
        #ax.set_xlabel("Opacity")
        ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.tick_params(axis='y', colors='white')
        ax.set_ylabel(label, color='white')

        
    def sample_colormap(self, v, w, alpha=None, colormap="gist_stern", col_bounds=None):
        r"""Add a Gaussian based on an existing colormap.

        Constructing pleasing Gaussians in a transfer function can pose some
        challenges, so this function will add a single Gaussian whose colors
        are taken from a colormap scaled between the bounds of the transfer
        function.  As with `TransferFunction.add_gaussian`, the value is
        calculated as :math:`f(x) = h \exp{-(x-x_0)^2 / w}` but with the height
        for each color calculated from the colormap.

        Parameters
        ----------
        v : float
            The value at which the Gaussian is to be added.
        w : float
            The relative width (:math:`w` in the above equation.)
        alpha : float, optional
            The alpha value height for the Gaussian
        colormap : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        col_bounds: array_like [min, max], optional
            Limits the values over which the colormap spans to these
            values.  Useful for sampling an entire colormap over a
            range smaller than the transfer function bounds.

        See Also
        --------
        ColorTransferFunction.add_layers : Many-at-a-time adder

        Examples
        --------

        >>> tf = ColorTransferFunction( (-10.0, -5.0) )
        >>> tf.sample_colormap(-7.0, 0.01, 'algae')
        """
        if col_bounds is None:
            rel = (v - self.x_bounds[0])/(self.x_bounds[1] - self.x_bounds[0])
        else:
            rel = (v - col_bounds[0])/(col_bounds[1] - col_bounds[0])
        cmap = get_cmap(colormap)
        r,g,b,a = cmap(rel)
        if alpha is None: alpha = a
        r *= alpha/3.
        g *= alpha/3.
        b *= alpha/3. 
        self.add_gaussian(v, w, [r, g, b, alpha])
        mylog.debug("Adding gaussian at %s with width %s and colors %s" % (
                v, w, (r,g,b,alpha)))

    def map_to_colormap(self, mi, ma, scale=1.0, colormap="gist_stern",
            scale_func=None):
        rel0 = int(self.nbins*(mi - self.x_bounds[0])/(self.x_bounds[1] -
            self.x_bounds[0]))
        rel1 = int(self.nbins*(ma - self.x_bounds[0])/(self.x_bounds[1] -
            self.x_bounds[0]))
        tomap = na.linspace(0.,1.,num=rel1-rel0)
        cmap = get_cmap(colormap)
        cc = cmap(tomap)*scale
        if scale_func is None:
            scale_mult = 1.0
        else:
            scale_mult = scale_func(tomap,0.0,1.0)
        print scale_mult 
        self.red.y[rel0:rel1]  = cc[:,0]*scale_mult
        self.green.y[rel0:rel1]= cc[:,1]*scale_mult
        self.blue.y[rel0:rel1] = cc[:,2]*scale_mult
        self.alpha.y[rel0:rel1]= cc[:,3]*scale_mult

    def add_layers(self, N, w=None, mi=None, ma=None, alpha = None,
                   colormap="gist_stern", col_bounds = None):
        r"""Add a set of Gaussians based on an existing colormap.

        Constructing pleasing Gaussians in a transfer function can pose some
        challenges, so this function will add several evenly-spaced Gaussians
        whose colors are taken from a colormap scaled between the bounds of the
        transfer function.   For each Gaussian to be added,
        `ColorTransferFunction.sample_colormap` is called.

        Parameters
        ----------
        N : int
            How many Gaussians to add
        w : float
            The relative width of each Gaussian.  If not included, it is
            calculated as 0.001 * (max_val - min_val) / N
        mi : float, optional
            If only a subset of the data range is to have the Gaussians added,
            this is the minimum for that subset
        ma : float, optional
            If only a subset of the data range is to have the Gaussians added,
            this is the maximum for that subset
        alpha : list of floats, optional
            The alpha value height for each Gaussian.  If not supplied, it is
            calculated as the logspace between -2.0 and 0.0.
        colormap : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        col_bounds: array_like [min, max], optional
            Limits the values over which the colormap spans to these
            values.  Useful for sampling an entire colormap over a
            range smaller than the transfer function bounds.

        See Also
        --------
        ColorTransferFunction.sample_colormap : Single Gaussian adder

        Examples
        --------

        >>> tf = ColorTransferFunction( (-10.0, -5.0) )
        >>> tf.add_layers(8)
        """
        if col_bounds is None:
            dist = (self.x_bounds[1] - self.x_bounds[0])
            if mi is None: mi = self.x_bounds[0] + dist/(10.0*N)
            if ma is None: ma = self.x_bounds[1] - dist/(10.0*N)
        else:
            dist = (col_bounds[1] - col_bounds[0])
            if mi is None: mi = col_bounds[0] + dist/(10.0*N)
            if ma is None: ma = col_bounds[1] - dist/(10.0*N)
        if w is None: w = 0.001 * (ma-mi)/N
        if alpha is None: alpha = na.logspace(-3.0, 0.0, N)
        for v, a in zip(na.mgrid[mi:ma:N*1j], alpha):
            self.sample_colormap(v, w, a, colormap=colormap, col_bounds=col_bounds)

    def get_colormap_image(self, height, width):
        image = na.zeros((height, width, 3), dtype='uint8')
        hvals = na.mgrid[self.x_bounds[0]:self.x_bounds[1]:height * 1j]
        for i,f in enumerate(self.funcs[:3]):
            vals = na.interp(hvals, f.x, f.y)
            image[:,:,i] = (vals[:,None] * 255).astype('uint8')
        image = image[::-1,:,:]
        return image
            
    def clear(self):
        for f in self.funcs:
            f.clear()


class ProjectionTransferFunction(MultiVariateTransferFunction):
    def __init__(self, x_bounds = (-1e60, 1e60), n_fields = 1):
        r"""A transfer function that defines a simple projection.

        To generate an interpolated, off-axis projection through a dataset,
        this transfer function should be used.  It will create a very simple
        table that merely sums along each ray.  Note that the end product will
        need to be scaled by the total width through which the rays were cast,
        a piece of information inacessible to the transfer function.

        Parameters
        ----------
        x_bounds : tuple of floats, optional
            If any of your values lie outside this range, they will be
            truncated.
        n_fields : int, optional
            How many fields we're going to project and pass through

        Notes
        -----
        When you use this transfer function, you may need to explicitly disable
        logging of fields.

        """
        if n_fields > 3:
            raise NotImplementedError
        MultiVariateTransferFunction.__init__(self)
        self.x_bounds = x_bounds
        self.nbins = 2
        self.linear_mapping = TransferFunction(x_bounds, 2)
        self.linear_mapping.pass_through = 1
        self.link_channels(0, [0,1,2]) # same emission for all rgb, default
        for i in range(n_fields):
            self.add_field_table(self.linear_mapping, i)
            self.link_channels(i, i)
        self.link_channels(n_fields, [3,4,5]) # this will remove absorption

class PlanckTransferFunction(MultiVariateTransferFunction):
    def __init__(self, T_bounds, rho_bounds, nbins=256,
                 red='R', green='V', blue='B',
                 anorm = 1e6):
        """
        This sets up a planck function for multivariate emission and
        absorption.  We assume that the emission is black body, which is then
        convolved with appropriate Johnson filters for *red*, *green* and
        *blue*.  *T_bounds* and *rho_bounds* define the limits of tabulated
        emission and absorption functions.  *anorm* is a "fudge factor" that
        defines the somewhat arbitrary normalization to the scattering
        approximation: because everything is done largely unit-free, and is
        really not terribly accurate anyway, feel free to adjust this to change
        the relative amount of reddenning.  Maybe in some future version this
        will be unitful.
        """
        MultiVariateTransferFunction.__init__(self)
        mscat = -1
        from UBVRI import johnson_filters
        for i, f in enumerate([red, green, blue]):
            jf = johnson_filters[f]
            tf = TransferFunction(T_bounds)
            tf.add_filtered_planck(jf['wavelen'], jf['trans'])
            self.add_field_table(tf, 0, 1)
            self.link_channels(i, i) # 0 => 0, 1 => 1, 2 => 2
            mscat = max(mscat, jf["Lchar"]**-4)

        for i, f in enumerate([red, green, blue]):
            # Now we set up the scattering
            scat = (johnson_filters[f]["Lchar"]**-4 / mscat)*anorm
            tf = TransferFunction(rho_bounds)
            mylog.debug("Adding: %s with relative scattering %s" % (f, scat))
            tf.y *= 0.0; tf.y += scat
            self.add_field_table(tf, 1, weight_field_id = 1)
            self.link_channels(i+3, i+3)

        self._normalize()

    def _normalize(self):
        fmax  = na.array([f.y for f in self.tables[:3]])
        normal = fmax.max(axis=0)
        for f in self.tables[:3]:
            f.y = f.y/normal

if __name__ == "__main__":
    tf = ColorTransferFunction((-20, -5))
    tf.add_gaussian(-16.0, 0.4, [0.2, 0.3, 0.1])
    tf.add_gaussian(-14.0, 0.8, [0.4, 0.1, 0.2])
    tf.add_gaussian(-10.0, 1.0, [0.0, 0.0, 1.0])
    tf.plot("tf.png")
