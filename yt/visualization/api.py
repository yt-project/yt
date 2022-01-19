from .base_plot_types import get_multi_plot
from .color_maps import add_colormap, make_colormap, show_colormaps
from .fits_image import (
    FITSImageData,
    FITSOffAxisProjection,
    FITSOffAxisSlice,
    FITSParticleProjection,
    FITSProjection,
    FITSSlice,
)
from .fixed_resolution import FixedResolutionBuffer, ParticleImageBuffer
from .image_writer import (
    apply_colormap,
    map_to_colors,
    multi_image_composite,
    scale_image,
    splat_points,
    write_bitmap,
    write_image,
    write_projection,
)
from .line_plot import LineBuffer, LinePlot
from .particle_plots import ParticlePhasePlot, ParticlePlot, ParticleProjectionPlot
from .plot_modifications import PlotCallback, callback_registry
from .plot_window import (
    AxisAlignedProjectionPlot,
    AxisAlignedSlicePlot,
    OffAxisProjectionPlot,
    OffAxisSlicePlot,
    ProjectionPlot,
    SlicePlot,
    plot_2d,
)
from .profile_plotter import PhasePlot, ProfilePlot
from .streamlines import Streamlines
