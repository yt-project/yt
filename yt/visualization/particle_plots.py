"""
This is a simple mechanism for interfacing with Particle plots



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import __builtin__
import base64
import types

from functools import wraps
from itertools import izip
import matplotlib
import numpy as np
import cStringIO

from yt.utilities.exceptions import \
    YTNotInsideNotebook
from yt.utilities.logger import ytLogger as mylog
import _mpl_imports as mpl
from yt.funcs import \
    ensure_list, \
    get_image_suffix, \
    get_ipython_api_version
from yt.units.unit_object import Unit
from .plot_window import \
    AxisAlignedParticlePlot
from .profile_plotter import \
    ParticlePhasePlot


def ParticlePlot(ds, x_field, y_field, z_fields=None, color='b', *args, **kwargs):
    r"""
    A factory function for
    :class:`yt.visualization.plot_window.AxisAlignedParticlePlot`
    and :class:`yt.visualization.profile_plotter.ParticlePhasePlot` objects. This
    essentially allows for a single entry point to both types of particle plots,
    the distinction being determined by the fields passed in. 

    If the x_field and y_field combination corresponds to a valid, right-handed
    spatial plot, an 'AxisAlignedParticlePlot` will be returned. This plot object 
    can be updated using one of the many helper functions defined in PlotWindow.

    If the x_field and y_field combo do not correspond to a valid 
    'AxisAlignedParticlePlot`, then a `ParticlePhasePlot`. This object can be
    modified by its own set of  helper functions defined in PhasePlot.

    Parameters
    ----------

    ds : :class:`yt.data_objects.api.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    x_field : string
        This is the particle field that will be plotted on the x-axis.
    y_field : string
        This is the particle field that will be plotted on the y-axis.
    z_fields : string, list, or None.
        If None, particles will be splatted onto the plot, but no colormap
        will be used. The particle color will instead be determined by
        the 'color' argument. If str or list, the name of the field or fields
        to be displayed on the colorbar.
        Default: None.
    color : 'b', 'g', 'r', 'c', 'm', 'y', 'k', or 'w'
         One the matplotlib-recognized color strings.
         The color that will indicate the particle locations
         on the plot. This argument is ignored if z_fields is
         not None. Default is 'b'.

    Examples
    --------

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y',
    ...                     'particle_mass', width=(0.5, 0.5))
    >>> p.set_unit('particle_mass', 'Msun')
    >>> p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_velocity_z',
    ...                     color='g')

    """

    direction = 3
    # try potential axes for a AxisAlignedParticlePlot:
    for axis in [0, 1, 2]:
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        ax_field_template = 'particle_position_%s'
        xf = ax_field_template % ds.coordinates.axis_name[xax]
        yf = ax_field_template % ds.coordinates.axis_name[yax]
        if (x_field, y_field) == (xf, yf):
            direction = axis
            break

    if direction < 3:
        # Make an AxisAlignedSlicePlot
        return AxisAlignedParticlePlot(ds, direction, z_fields, color,
                                       *args, **kwargs)

    # Does not correspond to any valid PlotWindow-style plot,
    # use ParticlePhasePlot instead
    else:
        return ParticlePhasePlot(ds.all_data(), x_field, y_field,
                                 z_fields, color, *args, **kwargs)
