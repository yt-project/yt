"""
API for yt.visualization



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .color_maps import \
    add_cmap, \
    show_colormaps, \
    add_cmap, \
    make_colormap

from .particle_plots import \
    ParticleProjectionPlot, \
    ParticlePhasePlot, \
    ParticlePlot

from .fixed_resolution import \
    FixedResolutionBuffer, \
    ObliqueFixedResolutionBuffer, \
    ParticleImageBuffer

from .image_writer import \
    multi_image_composite, \
    write_bitmap, \
    write_image, \
    map_to_colors, \
    splat_points, \
    apply_colormap, \
    scale_image, \
    write_projection

from .plot_modifications import \
    PlotCallback, \
    callback_registry

from .streamlines import \
    Streamlines

from .plot_window import \
    SlicePlot, \
    AxisAlignedSlicePlot, \
    OffAxisSlicePlot, \
    ProjectionPlot, \
    OffAxisProjectionPlot

from .profile_plotter import \
    ProfilePlot, \
    PhasePlot

from .base_plot_types import \
    get_multi_plot
