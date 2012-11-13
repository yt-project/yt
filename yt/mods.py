"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

from __future__ import absolute_import

#
# ALL IMPORTS GO HERE
#

# First module imports
import sys, types, os, glob, cPickle, time
import numpy as na # For historical reasons
import numpy as np # For modern purposes
import numpy # In case anyone wishes to use it by name

# This next item will handle most of the actual startup procedures, but it will
# also attempt to parse the command line and set up the global state of various
# operations.

import yt.startup_tasks as __startup_tasks
unparsed_args = __startup_tasks.unparsed_args

from yt.funcs import *
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.performance_counters import yt_counters, time_function
from yt.config import ytcfg, ytcfgDefaults
import yt.utilities.physical_constants as physical_constants

from yt.utilities.logger import level as __level
if __level >= int(ytcfgDefaults["loglevel"]):
    # This won't get displayed.
    mylog.debug("Turning off NumPy error reporting")
    np.seterr(all = 'ignore')

from yt.data_objects.api import \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    data_object_registry, \
    derived_field, add_field, FieldInfo, \
    ValidateParameter, ValidateDataField, ValidateProperty, \
    ValidateSpatial, ValidateGridType, \
    TimeSeriesData, AnalysisTask, analysis_task, \
    ParticleTrajectoryCollection, ImageArray

from yt.data_objects.derived_quantities import \
    add_quantity, quantity_info

from yt.frontends.enzo.api import \
    EnzoStaticOutput, EnzoStaticOutputInMemory, \
    EnzoSimulation, EnzoFieldInfo, \
    add_enzo_field, add_enzo_1d_field, add_enzo_2d_field

from yt.frontends.castro.api import \
    CastroStaticOutput, CastroFieldInfo, add_castro_field

from yt.frontends.nyx.api import \
    NyxStaticOutput, NyxFieldInfo, add_nyx_field

from yt.frontends.orion.api import \
    OrionStaticOutput, OrionFieldInfo, add_orion_field

from yt.frontends.flash.api import \
    FLASHStaticOutput, FLASHFieldInfo, add_flash_field

from yt.frontends.tiger.api import \
    TigerStaticOutput, TigerFieldInfo, add_tiger_field

from yt.frontends.artio.api import \
    ARTIOStaticOutput, ARTIOFieldInfo, add_artio_field

from yt.frontends.ramses.api import \
    RAMSESStaticOutput, RAMSESFieldInfo, add_ramses_field

from yt.frontends.chombo.api import \
    ChomboStaticOutput, ChomboFieldInfo, add_chombo_field

from yt.frontends.gdf.api import \
    GDFStaticOutput, GDFFieldInfo, add_gdf_field

from yt.frontends.athena.api import \
    AthenaStaticOutput, AthenaFieldInfo, add_athena_field

#from yt.frontends.art.api import \
#    ARTStaticOutput, ARTFieldInfo, add_art_field

#from yt.frontends.maestro.api import \
#    MaestroStaticOutput, MaestroFieldInfo, add_maestro_field

from yt.analysis_modules.list_modules import \
    get_available_modules, amods
available_analysis_modules = get_available_modules()

# Import our analysis modules
from yt.analysis_modules.halo_finding.api import \
    HaloFinder

from yt.utilities.definitions import \
    axis_names, x_dict, y_dict, inv_axis_names

# Now individual component imports from the visualization API
from yt.visualization.api import \
    PlotCollection, PlotCollectionInteractive, \
    get_multi_plot, FixedResolutionBuffer, ObliqueFixedResolutionBuffer, \
    callback_registry, write_bitmap, write_image, annotate_image, \
    apply_colormap, scale_image, write_projection, write_fits, \
    SlicePlot, OffAxisSlicePlot, ProjectionPlot, OffAxisProjectionPlot

from yt.visualization.volume_rendering.api import \
    ColorTransferFunction, PlanckTransferFunction, ProjectionTransferFunction, \
    HomogenizedVolume, Camera, off_axis_projection, MosaicFisheyeCamera

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects

for name, cls in callback_registry.items():
    exec("%s = cls" % name)

from yt.convenience import \
    load, projload, simulation

# Import some helpful math utilities
from yt.utilities.math_utils import \
    ortho_find, quartiles


# We load plugins.  Keep in mind, this can be fairly dangerous -
# the primary purpose is to allow people to have a set of functions
# that get used every time that they don't have to *define* every time.
# This way, other command-line tools can be used very simply.
# Unfortunately, for now, I think the easiest and simplest way of doing
# this is also the most dangerous way.
if ytcfg.getboolean("yt","loadfieldplugins"):
    my_plugin_name = ytcfg.get("yt","pluginfilename")
    # We assume that it is with respect to the $HOME/.yt directory
    if os.path.isfile(my_plugin_name):
        _fn = my_plugin_name
    else:
        _fn = os.path.expanduser("~/.yt/%s" % my_plugin_name)
    if os.path.isfile(_fn):
        mylog.info("Loading plugins from %s", _fn)
        execfile(_fn)
