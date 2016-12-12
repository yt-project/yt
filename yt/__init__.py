"""
YT is a package written primarily in Python designed to make the task of
running Enzo easier.  It contains facilities for creating Enzo data (currently
in prototype form) as well as runnning Enzo simulations, simulating the actions
of Enzo on various existing data, and analyzing output from Enzo in a
wide-variety of methods.

An ever-growing collection of documentation is also available at
http://yt-project.org/doc/ . Additionally, there is a
project site at http://yt-project.org/ with recipes, a wiki, a variety of
ways of peering into the version control, and a bug-reporting system.

YT is divided into several packages.

frontends
---------

This is where interfaces to codes are created.  Within each subdirectory of
yt/frontends/ there must exist the following files, even if empty:

* data_structures.py, where subclasses of AMRGridPatch, Dataset and
  AMRHierarchy are defined.
* io.py, where a subclass of IOHandler is defined.
* misc.py, where any miscellaneous functions or classes are defined.
* definitions.py, where any definitions specific to the frontend are
  defined.  (i.e., header formats, etc.)

visualization
-------------

This is where all visualization modules are stored.  This includes plot
collections, the volume rendering interface, and pixelization frontends.

data_objects
------------

All objects that handle data, processed or unprocessed, not explicitly
defined as visualization are located in here.  This includes the base
classes for data regions, covering grids, time series, and so on.  This
also includes derived fields and derived quantities.

analysis_modules
----------------

This is where all mechanisms for processing data live.  This includes
things like clump finding, halo profiling, halo finding, and so on.  This
is something of a catchall, but it serves as a level of greater
abstraction that simply data selection and modification.

gui
---

This is where all GUI components go.  Typically this will be some small
tool used for one or two things, which contains a launching mechanism on
the command line.

utilities
---------

All broadly useful code that doesn't clearly fit in one of the other
categories goes here.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = "3.3.3"

# First module imports
import numpy as np # For modern purposes
import numpy # In case anyone wishes to use it by name

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

import yt.utilities.physical_constants as physical_constants
import yt.units as units
from yt.units.yt_array import \
    YTArray, \
    YTQuantity, \
    uconcatenate, \
    uintersect1d, \
    uunion1d, \
    loadtxt, \
    savetxt

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
    DatasetSeries, ImageArray, \
    particle_filter, add_particle_filter, \
    create_profile, Profile1D, Profile2D, Profile3D, \
    ParticleProfile

# For backwards compatibility
TimeSeriesData = deprecated_class(DatasetSeries)

from yt.frontends.api import _frontend_container
frontends = _frontend_container()

from yt.frontends.stream.api import \
    load_uniform_grid, load_amr_grids, \
    load_particles, load_hexahedral_mesh, load_octree, \
    hexahedral_connectivity, load_unstructured_mesh

from yt.frontends.ytdata.api import \
    save_as_dataset

# For backwards compatibility
GadgetDataset = frontends.gadget.GadgetDataset
GadgetStaticOutput = deprecated_class(GadgetDataset)
TipsyDataset = frontends.tipsy.TipsyDataset
TipsyStaticOutput = deprecated_class(TipsyDataset)

# Now individual component imports from the visualization API
from yt.visualization.api import \
    FixedResolutionBuffer, ObliqueFixedResolutionBuffer, \
    write_bitmap, write_image, \
    apply_colormap, scale_image, write_projection, \
    SlicePlot, AxisAlignedSlicePlot, OffAxisSlicePlot, \
    ProjectionPlot, OffAxisProjectionPlot, \
    show_colormaps, add_cmap, make_colormap, \
    ProfilePlot, PhasePlot, ParticlePhasePlot, \
    ParticleProjectionPlot, ParticleImageBuffer, ParticlePlot

from yt.visualization.volume_rendering.api import \
    volume_render, create_scene, ColorTransferFunction, TransferFunction, \
    off_axis_projection, interactive_render
import yt.visualization.volume_rendering.api as volume_rendering
#    TransferFunctionHelper, MultiVariateTransferFunction
#    off_axis_projection

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, enable_parallelism, communication_system

from yt.convenience import \
    load, simulation

from yt.testing import run_nose

# Import some helpful math utilities
from yt.utilities.math_utils import \
    ortho_find, quartiles, periodic_position

from yt.units.unit_systems import UnitSystem
from yt.units.unit_object import unit_system_registry

from yt.analysis_modules.list_modules import \
    amods
