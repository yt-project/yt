"""
yt is a toolkit for analyzing and visualizing volumetric data.

* Website: http://yt-project.org
* Documentation: http://yt-project.org/doc
* Data hub: http://hub.yt
* Contribute: http://github.com/yt-project/yt

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = "3.4.1"

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
    deprecated_class, \
    toggle_interactivity
from yt.utilities.logger import ytLogger as mylog

import yt.utilities.physical_constants as physical_constants
import yt.units as units
from yt.units.unit_object import define_unit
from yt.units.yt_array import \
    YTArray, \
    YTQuantity, \
    uconcatenate, \
    ucross, \
    uintersect1d, \
    uunion1d, \
    unorm, \
    udot, \
    uvstack, \
    uhstack, \
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
    derived_field, \
    add_xray_emissivity_field

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
    SlicePlot, AxisAlignedSlicePlot, OffAxisSlicePlot, LinePlot, \
    LineBuffer, ProjectionPlot, OffAxisProjectionPlot, \
    show_colormaps, add_cmap, make_colormap, \
    ProfilePlot, PhasePlot, ParticlePhasePlot, \
    ParticleProjectionPlot, ParticleImageBuffer, ParticlePlot, \
    FITSImageData, FITSSlice, FITSProjection, FITSOffAxisSlice, \
    FITSOffAxisProjection, plot_2d

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
