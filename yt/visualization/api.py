"""
API for yt.visualization

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: J.S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

from color_maps import \
    add_cmap

from plot_collection import \
    PlotCollection, \
    PlotCollectionInteractive, \
    concatenate_pdfs, \
    get_multi_plot

from fixed_resolution import \
    FixedResolutionBuffer, \
    ObliqueFixedResolutionBuffer

from image_writer import \
    multi_image_composite, \
    write_bitmap, \
    write_image, \
    map_to_colors, \
    splat_points, \
    annotate_image, \
    apply_colormap, \
    scale_image, \
    write_projection, \
    write_fits

from plot_modifications import \
    PlotCallback, \
    callback_registry

from easy_plots import \
    plot_type_registry

from streamlines import \
    Streamlines

from plot_window import \
    SlicePlot, \
    OffAxisSlicePlot, \
    ProjectionPlot, \
    OffAxisProjectionPlot
    

