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
