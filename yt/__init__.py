"""
yt is a toolkit for analyzing and visualizing volumetric data.

* Website: https://yt-project.org
* Documentation: https://yt-project.org/doc
* Data hub: https://girder.hub.yt
* Contribute: https://github.com/yt-project/yt

"""
import sys

if sys.version_info[0] < 3:
    raise Exception(
        "Python 2 no longer supported.  Please install Python 3 for use with yt."
    )

__version__ = "4.0.dev0"

# First module imports
import numpy  # In case anyone wishes to use it by name
import numpy as np  # In case anyone wishes to use it by name

import yt.units as units
import yt.utilities.physical_constants as physical_constants
from yt.data_objects.api import (DatasetSeries, ImageArray, ParticleProfile,
                                 Profile1D, Profile2D, Profile3D,
                                 add_particle_filter, create_profile,
                                 particle_filter)
from yt.fields.api import (DerivedField, FieldDetector, FieldInfoContainer,
                           ValidateDataField, ValidateGridType,
                           ValidateParameter, ValidateProperty,
                           ValidateSpatial, add_field,
                           add_xray_emissivity_field, derived_field,
                           field_plugins)
from yt.funcs import (deprecated_class, enable_plugins, get_memory_usage,
                      get_pbar, get_version_stack, get_yt_supp, get_yt_version,
                      insert_ipython, is_root, iterable, memory_checker,
                      only_on_root, parallel_profile, print_tb, rootonly,
                      toggle_interactivity)
from yt.units import (YTArray, YTQuantity, display_ytarray, loadtxt, savetxt,
                      uconcatenate, ucross, udot, uhstack, uintersect1d, unorm,
                      ustack, uunion1d, uvstack)
from yt.units.unit_object import define_unit
from yt.utilities.logger import ytLogger as mylog

# For backwards compatibility
TimeSeriesData = deprecated_class(DatasetSeries)

from yt.frontends.api import _frontend_container

frontends = _frontend_container()

from yt.frontends.stream.api import (hexahedral_connectivity, load_amr_grids,
                                     load_hexahedral_mesh, load_octree,
                                     load_particles, load_uniform_grid,
                                     load_unstructured_mesh)
from yt.frontends.ytdata.api import save_as_dataset

# For backwards compatibility
GadgetDataset = frontends.gadget.GadgetDataset
GadgetStaticOutput = deprecated_class(GadgetDataset)
TipsyDataset = frontends.tipsy.TipsyDataset
TipsyStaticOutput = deprecated_class(TipsyDataset)

import yt.visualization.volume_rendering.api as volume_rendering
from yt.convenience import load, simulation
from yt.testing import run_nose
from yt.units.unit_systems import UnitSystem, unit_system_registry
from yt.utilities.load_sample import load_sample
# Import some helpful math utilities
from yt.utilities.math_utils import ortho_find, periodic_position, quartiles
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system, enable_parallelism, parallel_objects)
# Now individual component imports from the visualization API
from yt.visualization.api import (AxisAlignedSlicePlot, FITSImageData,
                                  FITSOffAxisProjection, FITSOffAxisSlice,
                                  FITSProjection, FITSSlice,
                                  FixedResolutionBuffer, LineBuffer, LinePlot,
                                  ObliqueFixedResolutionBuffer,
                                  OffAxisProjectionPlot, OffAxisSlicePlot,
                                  ParticleImageBuffer, ParticlePhasePlot,
                                  ParticlePlot, ParticleProjectionPlot,
                                  PhasePlot, ProfilePlot, ProjectionPlot,
                                  SlicePlot, add_colormap, apply_colormap,
                                  make_colormap, plot_2d, scale_image,
                                  show_colormaps, write_bitmap, write_image,
                                  write_projection)
from yt.visualization.volume_rendering.api import (ColorTransferFunction,
                                                   TransferFunction,
                                                   create_scene,
                                                   interactive_render,
                                                   off_axis_projection,
                                                   volume_render)

#    TransferFunctionHelper, MultiVariateTransferFunction
#    off_axis_projection







_called_from_pytest = False
