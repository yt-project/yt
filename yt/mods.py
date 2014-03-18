"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import

#
# ALL IMPORTS GO HERE
#

# First module imports
import numpy as np # For modern purposes
import numpy # In case anyone wishes to use it by name

# This next item will handle most of the actual startup procedures, but it will
# also attempt to parse the command line and set up the global state of various
# operations.  The variable unparsed_args is not used internally but is
# provided as a convenience for users who wish to parse arguments in scripts.
# See http://lists.spacepope.org/pipermail/yt-dev-spacepope.org/2011-December/
#     001727.html
import yt.startup_tasks as __startup_tasks
unparsed_args = __startup_tasks.unparsed_args

from yt.funcs import \
    iterable, \
    get_memory_usage, \
    print_tb, \
    rootonly, \
    insert_ipython, \
    get_pbar, \
    only_on_root, \
    is_root, \
    get_version_stack, \
    get_yt_supp, \
    get_yt_version, \
    parallel_profile, \
    enable_plugins, \
    memory_checker, \
    deprecated_class
from yt.utilities.logger import ytLogger as mylog
from yt.config import ytcfg, ytcfg_defaults
import yt.utilities.physical_constants as physical_constants
import yt.units as units
from yt.units.yt_array import YTArray, YTQuantity

from yt.utilities.logger import level as __level
if __level >= int(ytcfg_defaults["loglevel"]):
    # This won't get displayed.
    mylog.debug("Turning off NumPy error reporting")
    np.seterr(all = 'ignore')

from yt.fields.api import \
    field_plugins, \
    DerivedField, \
    FieldDetector, \
    FieldInfoContainer, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    add_field, \
    derived_field

from yt.data_objects.api import \
    BinnedProfile1D, BinnedProfile2D, BinnedProfile3D, \
    DatasetSeries, \
    ImageArray, particle_filter, create_profile, \
    Profile1D, Profile2D, Profile3D

from yt.frontends.api import _frontend_container
frontends = _frontend_container()

from yt.frontends.stream.api import \
    load_uniform_grid, load_amr_grids, \
    load_particles, load_hexahedral_mesh, load_octree

# For backwards compatibility
GadgetDataset = frontends.sph.GadgetDataset
GadgetStaticOutput = deprecated_class(GadgetDataset)
TipsyDataset = frontends.sph.TipsyDataset
TipsyStaticOutput = deprecated_class(TipsyDataset)

# Now individual component imports from the visualization API
from yt.visualization.api import \
    PlotCollection, PlotCollectionInteractive, \
    get_multi_plot, FixedResolutionBuffer, ObliqueFixedResolutionBuffer, \
    write_bitmap, write_image, \
    apply_colormap, scale_image, write_projection, \
    SlicePlot, AxisAlignedSlicePlot, OffAxisSlicePlot, \
    ProjectionPlot, OffAxisProjectionPlot, \
    show_colormaps, ProfilePlot, PhasePlot

from yt.visualization.volume_rendering.api import \
    off_axis_projection

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects

from yt.convenience import \
    load, simulation

# Import some helpful math utilities
from yt.utilities.math_utils import \
    ortho_find, quartiles, periodic_position


# We load plugins.  Keep in mind, this can be fairly dangerous -
# the primary purpose is to allow people to have a set of functions
# that get used every time that they don't have to *define* every time.
# This way, other command-line tools can be used very simply.
# Unfortunately, for now, I think the easiest and simplest way of doing
# this is also the most dangerous way.
if ytcfg.getboolean("yt","loadfieldplugins"):
    enable_plugins()
