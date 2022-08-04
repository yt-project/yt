"""
yt is a toolkit for analyzing and visualizing volumetric data.

* Website: https://yt-project.org
* Documentation: https://yt-project.org/doc
* Data hub: https://girder.hub.yt
* Contribute: https://github.com/yt-project/yt

"""
from ._version import __version__, version_info  # isort: skip
import yt.units as units
import yt.utilities.physical_constants as physical_constants
from yt.data_objects.api import (
    DatasetSeries,
    ImageArray,
    ParticleProfile,
    Profile1D,
    Profile2D,
    Profile3D,
    add_particle_filter,
    create_profile,
    particle_filter,
)
from yt.fields.api import (
    DerivedField,
    FieldDetector,
    FieldInfoContainer,
    ValidateDataField,
    ValidateGridType,
    ValidateParameter,
    ValidateProperty,
    ValidateSpatial,
    add_field,
    add_xray_emissivity_field,
    derived_field,
    field_plugins,
)
from yt.frontends.api import _frontend_container
from yt.funcs import (
    enable_plugins,
    get_memory_usage,
    get_pbar,
    get_version_stack,
    get_yt_version,
    insert_ipython,
    is_root,
    is_sequence,
    memory_checker,
    only_on_root,
    parallel_profile,
    print_tb,
    rootonly,
    toggle_interactivity,
)
from yt.units import (
    YTArray,
    YTQuantity,
    display_ytarray,
    loadtxt,
    savetxt,
    uconcatenate,
    ucross,
    udot,
    uhstack,
    uintersect1d,
    unorm,
    ustack,
    uunion1d,
    uvstack,
)
from yt.units.unit_object import define_unit  # type: ignore
from yt.utilities.logger import set_log_level, ytLogger as mylog

frontends = _frontend_container()

import yt.visualization.volume_rendering.api as volume_rendering
from yt.frontends.stream.api import hexahedral_connectivity
from yt.frontends.ytdata.api import save_as_dataset
from yt.loaders import (
    load,
    load_amr_grids,
    load_archive,
    load_hdf5_file,
    load_hexahedral_mesh,
    load_octree,
    load_particles,
    load_sample,
    load_simulation,
    load_uniform_grid,
    load_unstructured_mesh,
)


def run_nose(*args, **kwargs):
    # we hide this function behind a closure so we
    # don't make pytest a hard dependency for end users
    # see https://github.com/yt-project/yt/issues/3771
    from yt.testing import run_nose

    return run_nose(*args, **kwargs)


from yt.config import _setup_postinit_configuration
from yt.units.unit_systems import UnitSystem, unit_system_registry  # type: ignore

# Import some helpful math utilities
from yt.utilities.math_utils import ortho_find, periodic_position, quartiles
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system,
    enable_parallelism,
    parallel_objects,
)

# Now individual component imports from the visualization API
from yt.visualization.api import (
    AxisAlignedProjectionPlot,
    AxisAlignedSlicePlot,
    FITSImageData,
    FITSOffAxisProjection,
    FITSOffAxisSlice,
    FITSParticleProjection,
    FITSProjection,
    FITSSlice,
    FixedResolutionBuffer,
    LineBuffer,
    LinePlot,
    OffAxisProjectionPlot,
    OffAxisSlicePlot,
    ParticleImageBuffer,
    ParticlePhasePlot,
    ParticlePlot,
    ParticleProjectionPlot,
    PhasePlot,
    ProfilePlot,
    ProjectionPlot,
    SlicePlot,
    add_colormap,
    apply_colormap,
    make_colormap,
    plot_2d,
    scale_image,
    show_colormaps,
    write_bitmap,
    write_image,
    write_projection,
)
from yt.visualization.volume_rendering.api import (
    ColorTransferFunction,
    TransferFunction,
    create_scene,
    off_axis_projection,
    volume_render,
)

#    TransferFunctionHelper, MultiVariateTransferFunction
#    off_axis_projection


# run configuration callbacks
_setup_postinit_configuration()
del _setup_postinit_configuration
