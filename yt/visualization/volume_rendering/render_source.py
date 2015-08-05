"""
RenderSource Class

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog, ensure_numpy_array
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from .transfer_function_helper import TransferFunctionHelper
from .transfer_functions import TransferFunction, \
    ProjectionTransferFunction, ColorTransferFunction
from .utils import new_volume_render_sampler, data_source_or_all, \
    get_corners, new_projection_sampler, new_mesh_sampler
from yt.visualization.image_writer import apply_colormap
from .zbuffer_array import ZBuffer
from yt.utilities.lib.misc_utilities import \
    zlines, zpoints

from yt.utilities.on_demand_imports import NotAModule
try:
    from yt.utilities.lib import mesh_traversal
except ImportError:
    mesh_traversal = NotAModule("pyembree")
try:
    from yt.utilities.lib import mesh_construction
except ImportError:
    mesh_construction = NotAModule("pyembree")


class RenderSource(ParallelAnalysisInterface):

    """

    Base Class for Render Sources. Will be inherited for volumes,
    streamlines, etc.

    """

    def __init__(self):
        super(RenderSource, self).__init__()
        self.opaque = False
        self.zbuffer = None

    def render(self, camera, zbuffer=None):
        pass

    def _validate(self):
        pass


class OpaqueSource(RenderSource):
    """

    A base class for opaque render sources. Will be inherited from
    for LineSources, BoxSources, etc.

    """
    def __init__(self):
        super(OpaqueSource, self).__init__()
        self.opaque = True

    def set_zbuffer(self, zbuffer):
        self.zbuffer = zbuffer

    def render(self, camera, zbuffer=None):
        # This is definitely wrong for now
        if zbuffer is not None and self.zbuffer is not None:
            zbuffer.rgba = self.zbuffer.rgba
            zbuffer.z = self.zbuffer.z
            self.zbuffer = zbuffer
        return self.zbuffer


class VolumeSource(RenderSource):

    """

    A VolumeSource is a class for rendering data from
    an arbitrary volumetric data source, e.g. a sphere,
    cylinder, or the entire computational domain.


    """
    _image = None
    data_source = None

    def __init__(self, data_source, field, auto=True):
        r"""Initialize a new volumetric source for rendering.

        A :class:`VolumeSource` provides the framework to decompose an arbitrary
        yt data source into bricks that can be traversed and volume rendered.

        Parameters
        ----------
        data_source: :class:`AMR3DData` or :class:`Dataset`, optional
            This is the source to be rendered, which can be any arbitrary yt
            data object or dataset.
        fields : string
            The name of the field(s) to be rendered.
        auto: bool, optional
            If True, will build a default AMRKDTree and transfer function based
            on the data.

        Examples
        --------
        >>> source = RenderSource(ds, 'density')

        """
        super(VolumeSource, self).__init__()
        self.data_source = data_source_or_all(data_source)
        field = self.data_source._determine_fields(field)[0]
        self.field = field
        self.volume = None
        self.current_image = None
        self.double_check = False
        self.num_threads = 0
        self.num_samples = 10
        self.sampler_type = 'volume-render'

        # Error checking
        assert(self.field is not None)
        assert(self.data_source is not None)

        # In the future these will merge
        self.transfer_function = None
        self.tfh = None
        if auto:
            self.build_defaults()

    def build_defaults(self):
        self.build_default_volume()
        self.build_default_transfer_function()

    def set_transfer_function(self, transfer_function):
        """
        Set transfer function for this source
        """
        if not isinstance(transfer_function,
                          (TransferFunction, ColorTransferFunction,
                           ProjectionTransferFunction)):
            raise RuntimeError("transfer_function not of correct type")
        if isinstance(transfer_function, ProjectionTransferFunction):
            self.sampler_type = 'projection'

        self.transfer_function = transfer_function
        return self

    def _validate(self):
        """Make sure that all dependencies have been met"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.volume is None:
            raise RuntimeError("Volume not initialized")

        if self.transfer_function is None:
            raise RuntimeError("Transfer Function not Supplied")

    def build_default_transfer_function(self):
        self.tfh = \
            TransferFunctionHelper(self.data_source.pf)
        self.tfh.set_field(self.field)
        self.tfh.build_transfer_function()
        self.tfh.setup_default()
        self.transfer_function = self.tfh.tf

    def build_default_volume(self):
        self.volume = AMRKDTree(self.data_source.pf,
                                data_source=self.data_source)
        log_fields = [self.data_source.pf.field_info[self.field].take_log]
        mylog.debug('Log Fields:' + str(log_fields))
        self.volume.set_fields([self.field], log_fields, True)

    def set_volume(self, volume):
        assert(isinstance(volume, AMRKDTree))
        del self.volume
        self.volume = volume

    def set_field(self, field, no_ghost=True):
        field = self.data_source._determine_fields(field)[0]
        log_field = self.data_source.pf.field_info[field].take_log
        self.volume.set_fields(field, [log_field], no_ghost)
        self.field = field

    def set_fields(self, fields, no_ghost=True):
        fields = self.data_source._determine_fields(fields)
        log_fields = [self.data_source.ds.field_info[f].take_log
                      for f in fields]
        self.volume.set_fields(fields, log_fields, no_ghost)
        self.field = fields

    def set_sampler(self, camera):
        """docstring for add_sampler"""
        if self.sampler_type == 'volume-render':
            sampler = new_volume_render_sampler(camera, self)
        elif self.sampler_type == 'projection':
            sampler = new_projection_sampler(camera, self)
        else:
            NotImplementedError("%s not implemented yet" % self.sampler_type)
        self.sampler = sampler
        assert(self.sampler is not None)

    def render(self, camera, zbuffer=None):
        self.zbuffer = zbuffer
        self.set_sampler(camera)
        assert (self.sampler is not None)

        mylog.debug("Casting rays")
        total_cells = 0
        if self.double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(camera.lens.viewpoint):
            mylog.debug("Using sampler %s" % self.sampler)
            self.sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
        mylog.debug("Done casting rays")

        self.current_image = self.finalize_image(camera, self.sampler.aimage)
        if zbuffer is None:
            self.zbuffer = ZBuffer(self.current_image,
                                   np.inf*np.ones(self.current_image.shape[:2]))
        return self.current_image

    def finalize_image(self, camera, image):
        image = self.volume.reduce_tree_images(image,
                                               camera.lens.viewpoint)
        image.shape = camera.resolution[0], camera.resolution[1], 4
        if self.transfer_function.grey_opacity is False:
            image[:, :, 3] = 1.0
        return image

    def __repr__(self):
        disp = "<Volume Source>:%s " % str(self.data_source)
        disp += "transfer_function:%s" % str(self.transfer_function)
        return disp


class MeshSource(RenderSource):

    """docstring for MeshSource"""
    _image = None
    data_source = None

    def __init__(self, data_source, field):
        r"""Initialize a new unstructured source for rendering.

        A :class:`MeshSource` provides the framework to volume render
        unstructured mesh data.

        Parameters
        ----------
        data_source: :class:`AMR3DData` or :class:`Dataset`, optional
            This is the source to be rendered, which can be any arbitrary yt
            data object or dataset.
        fields : string
            The name of the field to be rendered.

        Examples
        --------
        >>> source = MeshSource(ds, ('all', 'convected'))

        """
        super(MeshSource, self).__init__()
        self.data_source = data_source_or_all(data_source)
        field = self.data_source._determine_fields(field)[0]
        self.field = field
        self.mesh = None
        self.current_image = None

        # Error checking
        assert(self.field is not None)
        assert(self.data_source is not None)

        self.scene = mesh_traversal.YTEmbreeScene()

        self.build_mesh()

    def _validate(self):
        """Make sure that all dependencies have been met"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

        if self.mesh is None:
            raise RuntimeError("Mesh not initialized")

    def build_mesh(self):

        field_data = self.data_source[self.field]
        vertices = self.data_source.ds.index.meshes[0].connectivity_coords

        # convert the indices to zero-based indexing
        indices = self.data_source.ds.index.meshes[0].connectivity_indices - 1

        self.mesh = mesh_construction.ElementMesh(self.scene,
                                                  vertices,
                                                  indices,
                                                  field_data.d)

    def render(self, camera):

        self.sampler = new_mesh_sampler(camera, self)

        mylog.debug("Casting rays")
        self.sampler(self.scene)
        mylog.debug("Done casting rays")

        self.current_image = self.sampler.aimage

        return self.current_image

    def __repr__(self):
        disp = "<Mesh Source>:%s " % str(self.data_source)
        return disp


class PointSource(OpaqueSource):

    _image = None
    data_source = None

    def __init__(self, positions, colors=None, color_stride=1):
        r"""A rendering source of opaque points in the scene.

        This class provides a mechanism for adding points to a scene; these
        points will be opaque, and can also be colored.

        Parameters
        ----------
        positions: array, shape (N, 3)
            These positions, in data-space coordinates, are the points to be
            added to the scene.
        colors : array, shape (N, 4), optional
            The colors of the points, including an alpha channel, in floating
            point running from 0..1.
        color_stride : int, optional
            The stride with which to access the colors when putting them on the
            scene.

        Examples
        --------
        >>> source = PointSource(particle_positions)

        """
        self.positions = positions
        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.ones((len(positions), 4))
            colors[:, 3] = 1.
        self.colors = colors
        self.color_stride = color_stride

    def render(self, camera, zbuffer=None):
        vertices = self.positions
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype='float64')
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # DRAW SOME POINTS
        camera.lens.setup_box_properties(camera)
        px, py, dz = camera.lens.project_to_plane(camera, vertices)
        zpoints(empty, z, px.d, py.d, dz.d, self.colors, self.color_stride)

        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Points Source>"
        return disp


class LineSource(OpaqueSource):

    _image = None
    data_source = None

    def __init__(self, positions, colors=None, color_stride=1):
        r"""A render source for a sequence of opaque line segments.

        This class provides a mechanism for adding lines to a scene; these
        points will be opaque, and can also be colored.

        Parameters
        ----------
        positions: array, shape (N, 2, 3)
            These positions, in data-space coordinates, are the starting and
            stopping points for each pair of lines. For example,
            positions[0][0] and positions[0][1] would give the (x, y, z)
            coordinates of the beginning and end points of the first line,
            respectively.
        colors : array, shape (N, 4), optional
            The colors of the points, including an alpha channel, in floating
            point running from 0..1.  Note that they correspond to the line
            segment succeeding each point; this means that strictly speaking
            they need only be (N-1) in length.
        color_stride : int, optional
            The stride with which to access the colors when putting them on the
            scene.

        Examples
        --------
        >>> source = LineSource(np.random.random((10, 3)))

        """

        super(LineSource, self).__init__()

        assert(positions.shape[1] == 2)
        assert(positions.shape[2] == 3)

        # convert the positions to the shape expected by zlines, below
        N = positions.shape[0]
        self.positions = positions.reshape((2*N, 3))

        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.ones((len(positions), 4))
            colors[:, 3] = 1.
        self.colors = colors
        self.color_stride = color_stride

    def render(self, camera, zbuffer=None):
        vertices = self.positions
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype='float64')
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # DRAW SOME LINES
        camera.lens.setup_box_properties(camera)
        px, py, dz = camera.lens.project_to_plane(camera, vertices)
        zlines(empty, z, px.d, py.d, dz.d, self.colors, self.color_stride)

        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Line Source>"
        return disp


class BoxSource(LineSource):
    def __init__(self, left_edge, right_edge, color=None):
        r"""A render source for a box drawn with line segments.

        This render source will draw a box, with transparent faces, in data
        space coordinates.  This is useful for annotations.

        Parameters
        ----------
        left_edge: array-like, shape (3,), float
            The left edge coordinates of the box.
        right_edge : array-like, shape (3,), float
            The right edge coordinates of the box.
        color : array-like, shape (4,), float, optional
            The colors (including alpha) to use for the lines.

        Examples
        --------
        >>> source = BoxSource(grid_obj.LeftEdge, grid_obj.RightEdge)

        """
        if color is None:
            color = np.array([1.0, 1.0, 1.0, 1.0])
        color = ensure_numpy_array(color)
        color.shape = (1, 4)
        corners = get_corners(left_edge.copy(), right_edge.copy())
        order = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]
        vertices = np.empty([24, 3])
        for i in range(3):
            vertices[:, i] = corners[order, i, ...].ravel(order='F')
        vertices = vertices.reshape((12, 2, 3))

        super(BoxSource, self).__init__(vertices, color, color_stride=24)


class GridSource(LineSource):
    def __init__(self, data_source, alpha=0.3, cmap='alage',
                 min_level=None, max_level=None):
        r"""A render source for drawing grids in a scene.

        This render source will draw blocks that are within a given data
        source, by default coloring them by their level of resolution.

        Parameters
        ----------
        data_source: :class:`~yt.data_objects.api.DataContainer`
            The data container that will be used to identify grids to draw.
        alpha : float
            The opacity of the grids to draw.
        cmap : color map name
            The color map to use to map resolution levels to color.
        min_level : int, optional
            Minimum level to draw
        max_level : int, optional
            Maximum level to draw

        Examples
        --------
        >>> dd = ds.sphere("c", (0.1, "unitary"))
        >>> source = GridSource(dd, alpha=1.0)

        """
        data_source = data_source_or_all(data_source)
        corners = []
        levels = []
        for block, mask in data_source.blocks:
            block_corners = np.array([
                [block.LeftEdge[0], block.LeftEdge[1], block.LeftEdge[2]],
                [block.RightEdge[0], block.LeftEdge[1], block.LeftEdge[2]],
                [block.RightEdge[0], block.RightEdge[1], block.LeftEdge[2]],
                [block.LeftEdge[0], block.RightEdge[1], block.LeftEdge[2]],
                [block.LeftEdge[0], block.LeftEdge[1], block.RightEdge[2]],
                [block.RightEdge[0], block.LeftEdge[1], block.RightEdge[2]],
                [block.RightEdge[0], block.RightEdge[1], block.RightEdge[2]],
                [block.LeftEdge[0], block.RightEdge[1], block.RightEdge[2]],
            ], dtype='float64')
            corners.append(block_corners)
            levels.append(block.Level)
        corners = np.dstack(corners)
        levels = np.array(levels)

        if max_level is not None:
            subset = levels <= max_level
            levels = levels[subset]
            corners = corners[:, :, subset]
        if min_level is not None:
            subset = levels >= min_level
            levels = levels[subset]
            corners = corners[:, :, subset]

        colors = apply_colormap(
            levels*1.0,
            color_bounds=[0, data_source.ds.index.max_level],
            cmap_name=cmap)[0, :, :]*alpha/255.
        colors[:, 3] = alpha

        order = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]

        vertices = np.empty([corners.shape[2]*2*12, 3])
        for i in range(3):
            vertices[:, i] = corners[order, i, ...].ravel(order='F')
        vertices = vertices.reshape((corners.shape[2]*12, 2, 3))

        super(GridSource, self).__init__(vertices, colors, color_stride=24)


class CoordinateVectorSource(OpaqueSource):
    def __init__(self, colors=None, alpha=1.0):
        r"""Draw coordinate vectors on the scene.

        This will draw a set of coordinate vectors on the camera image.  They
        will appear in the lower right of the image.

        Parameters
        ----------
        colors: array-like, shape (3,4), optional
            The x, y, z RGBA values to use to draw the vectors.
        alpha : float, optional
            The opacity of the vectors.

        Examples
        --------
        >>> source = CoordinateVectorSource()

        """
        super(CoordinateVectorSource, self).__init__()
        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.zeros((3, 4))
            colors[0, 0] = alpha  # x is red
            colors[1, 1] = alpha  # y is green
            colors[2, 2] = alpha  # z is blue
            colors[:, 3] = alpha
        self.colors = colors
        self.color_stride = 2

    def render(self, camera, zbuffer=None):
        camera.lens.setup_box_properties(camera)
        center = camera.focus
        # Get positions at the focus
        positions = np.zeros([6, 3])
        positions[:] = center

        # Create vectors in the x,y,z directions
        for i in range(3):
            positions[2*i+1, i] += camera.width.d[i] / 16.0

        # Project to the image plane
        px, py, dz = camera.lens.project_to_plane(camera, positions)
        dpx = px[1::2] - px[::2]
        dpy = py[1::2] - py[::2]

        # Set the center of the coordinates to be in the lower left of the image
        lpx = camera.resolution[0] / 8
        lpy = camera.resolution[1] - camera.resolution[1] / 8  # Upside-downsies

        # Offset the pixels according to the projections above
        px[::2] = lpx
        px[1::2] = lpx + dpx
        py[::2] = lpy
        py[1::2] = lpy + dpy
        dz[:] = 0.0

        # Create a zbuffer if needed
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype='float64')
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # Draw the vectors
        zlines(empty, z, px.d, py.d, dz.d, self.colors, self.color_stride)

        # Set the new zbuffer
        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Coordinates Source>"
        return disp
