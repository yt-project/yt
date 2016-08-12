"""
Import the components of the volume rendering extension



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.extern.six.moves import builtins
import numpy as np

from yt.config import \
    ytcfg
from yt.funcs import \
    iterable, mylog, get_pbar, \
    get_num_threads, ensure_numpy_array
from yt.units.yt_array import YTArray
from yt.utilities.exceptions import YTNotInsideNotebook
from copy import deepcopy

from .transfer_functions import ProjectionTransferFunction

from yt.utilities.lib.grid_traversal import \
    pixelize_healpix, arr_fisheye_vectors, arr_pix2vec_nest
from yt.utilities.lib.partitioned_grid import \
    PartitionedGrid
from yt.utilities.lib.image_samplers import \
    ProjectionSampler, VolumeRenderSampler, \
    LightSourceRenderSampler, InterpolatedProjectionSampler
from yt.utilities.lib.misc_utilities import \
    lines

from yt.utilities.math_utils import get_rotation_matrix
from yt.utilities.orientation import Orientation
from yt.data_objects.api import ImageArray
from yt.visualization.image_writer import write_bitmap, write_image, apply_colormap
from yt.data_objects.data_containers import data_object_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_objects
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.visualization.volume_rendering.blenders import enhance_rgba

def get_corners(le, re):
    return np.array([
      [le[0], le[1], le[2]],
      [re[0], le[1], le[2]],
      [re[0], re[1], le[2]],
      [le[0], re[1], le[2]],
      [le[0], le[1], re[2]],
      [re[0], le[1], re[2]],
      [re[0], re[1], re[2]],
      [le[0], re[1], re[2]],
    ], dtype='float64')

class Camera(ParallelAnalysisInterface):
    r"""A viewpoint into a volume, for volume rendering.

    The camera represents the eye of an observer, which will be used to
    generate ray-cast volume renderings of the domain.

    Parameters
    ----------
    center : array_like
        The current "center" of the view port -- the focal point for the
        camera.
    normal_vector : array_like
        The vector between the camera position and the center.
    width : float or list of floats
        The current width of the image.  If a single float, the volume is
        cubical, but if not, it is left/right, top/bottom, front/back.
    resolution : int or list of ints
        The number of pixels in each direction.
    transfer_function : `yt.visualization.volume_rendering.TransferFunction`
        The transfer function used to map values to colors in an image.  If
        not specified, defaults to a ProjectionTransferFunction.
    north_vector : array_like, optional
        The 'up' direction for the plane of rays.  If not specific, calculated
        automatically.
    steady_north : bool, optional
        Boolean to control whether to normalize the north_vector
        by subtracting off the dot product of it and the normal
        vector.  Makes it easier to do rotations along a single
        axis.  If north_vector is specified, is switched to
        True. Default: False
    volume : `yt.extensions.volume_rendering.AMRKDTree`, optional
        The volume to ray cast through.  Can be specified for finer-grained
        control, but otherwise will be automatically generated.
    fields : list of fields, optional
        This is the list of fields we want to volume render; defaults to
        Density.
    log_fields : list of bool, optional
        Whether we should take the log of the fields before supplying them to
        the volume rendering mechanism.
    sub_samples : int, optional
        The number of samples to take inside every cell per ray.
    ds : `~yt.data_objects.api.Dataset`
        For now, this is a require parameter!  But in the future it will become
        optional.  This is the dataset to volume render.
    use_kd: bool, optional
        Specifies whether or not to use a kd-Tree framework for
        the Homogenized Volume and ray-casting.  Default to True.
    max_level: int, optional
        Specifies the maximum level to be rendered.  Also
        specifies the maximum level used in the kd-Tree
        construction.  Defaults to None (all levels), and only
        applies if use_kd=True.
    no_ghost: bool, optional
        Optimization option.  If True, homogenized bricks will
        extrapolate out from grid instead of interpolating from
        ghost zones that have to first be calculated.  This can
        lead to large speed improvements, but at a loss of
        accuracy/smoothness in resulting image.  The effects are
        less notable when the transfer function is smooth and
        broad. Default: True
    data_source: data container, optional
        Optionally specify an arbitrary data source to the volume rendering.
        All cells not included in the data source will be ignored during ray
        casting. By default this will get set to ds.all_data().

    Examples
    --------

    >>> from yt.mods import *
    >>> import yt.visualization.volume_rendering.api as vr

    >>> ds = load('DD1701') # Load a dataset
    >>> c = [0.5]*3 # Center
    >>> L = [1.0,1.0,1.0] # Viewpoint
    >>> W = np.sqrt(3) # Width
    >>> N = 1024 # Pixels (1024^2)

    # Get density min, max
    >>> mi, ma = ds.all_data().quantities['Extrema']('Density')[0]
    >>> mi, ma = np.log10(mi), np.log10(ma)

    # Construct transfer function
    >>> tf = vr.ColorTransferFunction((mi-2, ma+2))
    # Sample transfer function with 5 gaussians.  Use new col_bounds keyword.
    >>> tf.add_layers(5,w=0.05, col_bounds = (mi+1,ma), colormap='spectral')
    
    # Create the camera object
    >>> cam = vr.Camera(c, L, W, (N,N), transfer_function=tf, ds=ds)
    
    # Ray cast, and save the image.
    >>> image = cam.snapshot(fn='my_rendering.png')

    """
    _sampler_object = VolumeRenderSampler
    _pylab = None
    _tf_figure = None
    _render_figure = None
    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function = None,
                 north_vector = None, steady_north=False,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, ds = None,
                 min_level=None, max_level=None, no_ghost=True,
                 data_source=None,
                 use_light=False):
        ParallelAnalysisInterface.__init__(self)
        if ds is not None: self.ds = ds
        if not iterable(resolution):
            resolution = (resolution, resolution)
        self.resolution = resolution
        self.sub_samples = sub_samples
        self.rotation_vector = north_vector
        if iterable(width) and len(width) > 1 and isinstance(width[1], str):
            width = self.ds.quan(width[0], input_units=width[1])
            # Now convert back to code length for subsequent manipulation
            width = width.in_units("code_length").value
        if not iterable(width):
            width = (width, width, width) # left/right, top/bottom, front/back 
        if not isinstance(width, YTArray):
            width = self.ds.arr(width, input_units="code_length")
        if not isinstance(center, YTArray):
            center = self.ds.arr(center, input_units="code_length")
        # Ensure that width and center are in the same units
        # Cf. https://bitbucket.org/yt_analysis/yt/issue/1080
        width.convert_to_units("code_length")
        center.convert_to_units("code_length")
        self.orienter = Orientation(normal_vector, north_vector=north_vector, steady_north=steady_north)
        if not steady_north:
            self.rotation_vector = self.orienter.unit_vectors[1]
        self._setup_box_properties(width, center, self.orienter.unit_vectors)
        if fields is None: fields = ["density"]
        self.fields = fields
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        self.log_fields = log_fields
        dd = self.ds.all_data()
        efields = dd._determine_fields(self.fields)
        if self.log_fields is None:
            self.log_fields = [self.ds._get_field_info(*f).take_log for f in efields]
        self.no_ghost = no_ghost
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if self.no_ghost:
            mylog.info('Warning: no_ghost is currently True (default). This may lead to artifacts at grid boundaries.')

        if data_source is None:
            data_source = self.ds.all_data()
        self.data_source = data_source

        if volume is None:
            volume = AMRKDTree(self.ds, min_level=min_level, 
                               max_level=max_level, data_source=self.data_source)
        self.volume = volume        

    def _setup_box_properties(self, width, center, unit_vectors):
        self.width = width
        self.center = center
        self.box_vectors = YTArray([unit_vectors[0]*width[0],
                                    unit_vectors[1]*width[1],
                                    unit_vectors[2]*width[2]])
        self.origin = center - 0.5*width.dot(YTArray(unit_vectors, ""))
        self.back_center =  center - 0.5*width[2]*unit_vectors[2]
        self.front_center = center + 0.5*width[2]*unit_vectors[2]         

    def update_view_from_matrix(self, mat):
        pass

    def project_to_plane(self, pos, res=None):
        if res is None: 
            res = self.resolution
        dx = np.dot(pos - self.origin, self.orienter.unit_vectors[1])
        dy = np.dot(pos - self.origin, self.orienter.unit_vectors[0])
        dz = np.dot(pos - self.center, self.orienter.unit_vectors[2])
        # Transpose into image coords.
        py = (res[0]*(dx/self.width[0])).astype('int')
        px = (res[1]*(dy/self.width[1])).astype('int')
        return px, py, dz

    def draw_grids(self, im, alpha=0.3, cmap=None, min_level=None, 
                   max_level=None):
        r"""Draws Grids on an existing volume rendering.

        By mapping grid level to a color, draws edges of grids on 
        a volume rendering using the camera orientation.

        Parameters
        ----------
        im: Numpy ndarray
            Existing image that has the same resolution as the Camera, 
            which will be painted by grid lines.
        alpha : float, optional
            The alpha value for the grids being drawn.  Used to control
            how bright the grid lines are with respect to the image.
            Default : 0.3
        cmap : string, optional
            Colormap to be used mapping grid levels to colors.
        min_level, max_level : int, optional
            Optional parameters to specify the min and max level grid boxes 
            to overplot on the image.  
        
        Returns
        -------
        None

        Examples
        --------
        >>> im = cam.snapshot() 
        >>> cam.add_grids(im)
        >>> write_bitmap(im, 'render_with_grids.png')

        """
        if cmap is None:
            cmap = ytcfg.get("yt", "default_colormap")
        region = self.data_source
        corners = []
        levels = []
        for block, mask in region.blocks:
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
            corners = corners[:,:,subset]
        if min_level is not None:
            subset = levels >= min_level
            levels = levels[subset]
            corners = corners[:,:,subset]
            
        colors = apply_colormap(levels*1.0,
                                color_bounds=[0,self.ds.index.max_level],
                                cmap_name=cmap)[0,:,:]*1.0/255.
        colors[:,3] = alpha

                
        order  = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]
        
        vertices = np.empty([corners.shape[2]*2*12,3])
        vertices = self.ds.arr(vertices, "code_length")
        for i in range(3):
            vertices[:,i] = corners[order,i,...].ravel(order='F')

        px, py, dz = self.project_to_plane(vertices, res=im.shape[:2])
        
        # Must normalize the image
        nim = im.rescale(inline=False)
        enhance_rgba(nim)
        nim.add_background_color('black', inline=True)

        # we flipped it in snapshot to get the orientation correct, so
        # flip the lines
        lines(nim.d, px.d, py.d, colors, 24, flip=1)

        return nim

    def draw_coordinate_vectors(self, im, length=0.05, thickness=1):
        r"""Draws three coordinate vectors in the corner of a rendering.

        Modifies an existing image to have three lines corresponding to the
        coordinate directions colored by {x,y,z} = {r,g,b}.  Currently only
        functional for plane-parallel volume rendering.

        Parameters
        ----------
        im: Numpy ndarray
            Existing image that has the same resolution as the Camera,
            which will be painted by grid lines.
        length: float, optional
            The length of the lines, as a fraction of the image size.
            Default : 0.05
        thickness : int, optional
            Thickness in pixels of the line to be drawn.

        Returns
        -------
        None

        Modifies
        --------
        im: The original image.

        Examples
        --------
        >>> im = cam.snapshot()
        >>> cam.draw_coordinate_vectors(im)
        >>> im.write_png('render_with_grids.png')

        """
        length_pixels = length * self.resolution[0]
        # Put the starting point in the lower left
        px0 = int(length * self.resolution[0])
        # CS coordinates!
        py0 = int((1.0-length) * self.resolution[1])

        alpha = im[:, :, 3].max()
        if alpha == 0.0:
            alpha = 1.0

        coord_vectors = [np.array([length_pixels, 0.0, 0.0]),
                         np.array([0.0, length_pixels, 0.0]),
                         np.array([0.0, 0.0, length_pixels])]
        colors = [np.array([1.0, 0.0, 0.0, alpha]),
                  np.array([0.0, 1.0, 0.0, alpha]),
                  np.array([0.0, 0.0, 1.0, alpha])]

        # we flipped it in snapshot to get the orientation correct, so
        # flip the lines
        for vec, color in zip(coord_vectors, colors):
            dx = int(np.dot(vec, self.orienter.unit_vectors[0]))
            dy = int(np.dot(vec, self.orienter.unit_vectors[1]))
            px = np.array([px0, px0+dx], dtype='int64')
            py = np.array([py0, py0+dy], dtype='int64')
            lines(im.d, px, py, np.array([color, color]), 1, thickness, flip=1)

    def draw_line(self, im, x0, x1, color=None):
        r"""Draws a line on an existing volume rendering.
        Given starting and ending positions x0 and x1, draws a line on 
        a volume rendering using the camera orientation.

        Parameters
        ----------
        im : ImageArray or 2D ndarray
            Existing image that has the same resolution as the Camera, 
            which will be painted by grid lines.
        x0 : YTArray or ndarray
            Starting coordinate.  If passed in as an ndarray,
            assumed to be in code units.
        x1 : YTArray or ndarray
            Ending coordinate, in simulation coordinates.  If passed in as
            an ndarray, assumed to be in code units.
        color : array like, optional
            Color of the line (r, g, b, a). Defaults to white. 

        Returns
        -------
        None

        Examples
        --------
        >>> im = cam.snapshot() 
        >>> cam.draw_line(im, np.array([0.1,0.2,0.3], np.array([0.5,0.6,0.7)))
        >>> write_bitmap(im, 'render_with_line.png')

        """
        if color is None:
            color = np.array([1.0,1.0,1.0,1.0])

        if not hasattr(x0, "units"):
            x0 = self.ds.arr(x0, "code_length")
        if not hasattr(x1, "units"):
            x1 = self.ds.arr(x1, "code_length")

        dx0 = ((x0-self.origin)*self.orienter.unit_vectors[1]).sum()
        dx1 = ((x1-self.origin)*self.orienter.unit_vectors[1]).sum()
        dy0 = ((x0-self.origin)*self.orienter.unit_vectors[0]).sum()
        dy1 = ((x1-self.origin)*self.orienter.unit_vectors[0]).sum()
        py0 = int(self.resolution[0]*(dx0/self.width[0]))
        py1 = int(self.resolution[0]*(dx1/self.width[0]))
        px0 = int(self.resolution[1]*(dy0/self.width[1]))
        px1 = int(self.resolution[1]*(dy1/self.width[1]))
        px = np.array([px0, px1], dtype="int64")
        py = np.array([py0, py1], dtype="int64")
        # we flipped it in snapshot to get the orientation correct, so
        # flip the lines
        lines(im.d, px, py, np.array([color,color]), flip=1)

    def draw_domain(self,im,alpha=0.3):
        r"""Draws domain edges on an existing volume rendering.

        Draws a white wireframe on the domain edges.

        Parameters
        ----------
        im: Numpy ndarray
            Existing image that has the same resolution as the Camera, 
            which will be painted by grid lines.
        alpha : float, optional
            The alpha value for the wireframe being drawn.  Used to control
            how bright the lines are with respect to the image.
            Default : 0.3
        
        Returns
        -------
        nim: Numpy ndarray
            A new image with the domain lines drawn

        Examples
        --------
        >>> im = cam.snapshot() 
        >>> nim = cam.draw_domain(im)
        >>> write_bitmap(nim, 'render_with_domain_boundary.png')

        """
        # Must normalize the image
        nim = im.rescale(inline=False)
        enhance_rgba(nim)
        nim.add_background_color('black', inline=True)
 
        self.draw_box(nim, self.ds.domain_left_edge, self.ds.domain_right_edge,
                        color=np.array([1.0,1.0,1.0,alpha]))
        return nim

    def draw_box(self, im, le, re, color=None):
        r"""Draws a box on an existing volume rendering.

        Draws a box defined by a left and right edge by modifying an
        existing volume rendering

        Parameters
        ----------
        im: Numpy ndarray
            Existing image that has the same resolution as the Camera, 
            which will be painted by grid lines.
        le: Numpy ndarray
            Left corner of the box 
        re : Numpy ndarray
            Right corner of the box 
        color : array like, optional
            Color of the box (r, g, b, a). Defaults to white. 
        
        Returns
        -------
        None

        Examples
        --------
        >>> im = cam.snapshot() 
        >>> cam.draw_box(im, np.array([0.1,0.2,0.3], np.array([0.5,0.6,0.7)))
        >>> write_bitmap(im, 'render_with_box.png')

        """

        if color is None:
            color = np.array([1.0,1.0,1.0,1.0]) 
        corners = get_corners(le,re)
        order  = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]
        
        vertices = np.empty([24,3])
        vertices = self.ds.arr(vertices, "code_length")
        for i in range(3):
            vertices[:,i] = corners[order,i,...].ravel(order='F')

        px, py, dz = self.project_to_plane(vertices, res=im.shape[:2])
       
        # we flipped it in snapshot to get the orientation correct, so
        # flip the lines
        lines(im.d, px.d.astype("int64"), py.d.astype("int64"), color.reshape(1,4), 24, flip=1)

    def look_at(self, new_center, north_vector = None):
        r"""Change the view direction based on a new focal point.

        This will recalculate all the necessary vectors and vector planes to orient
        the image plane so that it points at a new location.

        Parameters
        ----------
        new_center : array_like
            The new "center" of the view port -- the focal point for the
            camera.
        north_vector : array_like, optional
            The "up" direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        normal_vector = self.front_center - new_center
        self.orienter.switch_orientation(normal_vector=normal_vector,
                                         north_vector = north_vector)

    def switch_orientation(self, normal_vector=None, north_vector=None):
        r"""
        Change the view direction based on any of the orientation parameters.

        This will recalculate all the necessary vectors and vector planes
        related to an orientable object.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if north_vector is None:
            north_vector = self.north_vector
        if normal_vector is None:
            normal_vector = self.normal_vector
        self.orienter._setup_normalized_vectors(normal_vector, north_vector)

    def switch_view(self, normal_vector=None, width=None, center=None, north_vector=None):
        r"""Change the view based on any of the view parameters.

        This will recalculate the orientation and width based on any of
        normal_vector, width, center, and north_vector.

        Parameters
        ----------
        normal_vector: array_like, optional
            The new looking vector.
        width: float or array of floats, optional
            The new width.  Can be a single value W -> [W,W,W] or an
            array [W1, W2, W3] (left/right, top/bottom, front/back)
        center: array_like, optional
            Specifies the new center.
        north_vector : array_like, optional
            The 'up' direction for the plane of rays.  If not specific,
            calculated automatically.
        """
        if width is None:
            width = self.width
        if not iterable(width):
            width = (width, width, width) # left/right, tom/bottom, front/back 
        self.width = width
        if center is not None:
            self.center = center
        if north_vector is None:
            north_vector = self.orienter.north_vector
        if normal_vector is None:
            normal_vector = self.orienter.normal_vector
        self.switch_orientation(normal_vector = normal_vector,
                                         north_vector = north_vector)
        self._setup_box_properties(width, self.center, self.orienter.unit_vectors)
        
    def new_image(self):
        image = np.zeros((self.resolution[0], self.resolution[1], 4), dtype='float64', order='C')
        return image

    def get_sampler_args(self, image):
        rotp = np.concatenate([self.orienter.inv_mat.ravel('F'), self.back_center.ravel()])
        args = (np.atleast_3d(rotp), np.atleast_3d(self.box_vectors[2]),
                self.back_center,
                (-self.width[0]/2.0, self.width[0]/2.0,
                 -self.width[1]/2.0, self.width[1]/2.0),
                image, self.orienter.unit_vectors[0], self.orienter.unit_vectors[1],
                np.array(self.width, dtype='float64'), self.transfer_function, self.sub_samples)
        return args, {'lens_type': 'plane-parallel'}

    star_trees = None
    def get_sampler(self, args, kwargs):
        if self.star_trees is not None:
            kwargs = {'star_list': self.star_trees}
        if self.use_light:
            if self.light_dir is None:
                self.set_default_light_dir()
            temp_dir = np.empty(3,dtype='float64')
            temp_dir = self.light_dir[0] * self.orienter.unit_vectors[1] + \
                    self.light_dir[1] * self.orienter.unit_vectors[2] + \
                    self.light_dir[2] * self.orienter.unit_vectors[0]
            if self.light_rgba is None:
                self.set_default_light_rgba()
            sampler = LightSourceRenderSampler(*args, light_dir=temp_dir,
                    light_rgba=self.light_rgba, **kwargs)
        else:
            sampler = self._sampler_object(*args, **kwargs)
        return sampler

    def finalize_image(self, image):
        view_pos = self.front_center + self.orienter.unit_vectors[2] * 1.0e6 * self.width[2]
        image = self.volume.reduce_tree_images(image, view_pos)
        if self.transfer_function.grey_opacity is False:
            image[:,:,3]=1.0
        return image

    def _render(self, double_check, num_threads, image, sampler):
        ncells = sum(b.source_mask.size for b in self.volume.bricks)
        pbar = get_pbar("Ray casting", ncells)
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        view_pos = self.front_center + self.orienter.unit_vectors[2] * 1.0e6 * self.width[2]
        for brick in self.volume.traverse(view_pos):
            sampler(brick, num_threads=num_threads)
            total_cells += brick.source_mask.size
            pbar.update(total_cells)

        pbar.finish()
        image = sampler.aimage
        image = self.finalize_image(image)
        return image

    def show_tf(self):
        if self._pylab is None: 
            import pylab
            self._pylab = pylab
        if self._tf_figure is None:
            self._tf_figure = self._pylab.figure(2)
            self.transfer_function.show(ax=self._tf_figure.axes)
        self._pylab.draw()

    def annotate(self, ax, enhance=True, label_fmt=None):
        ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_visible(False)
        ax.get_yaxis().set_ticks([])
        cb = self._pylab.colorbar(ax.images[0], pad=0.0, fraction=0.05, drawedges=True, shrink=0.9)
        label = self.ds._get_field_info(self.fields[0]).get_label()
        if self.log_fields[0]:
            label = r'$\rm{log}\ $' + label
        self.transfer_function.vert_cbar(ax=cb.ax, label=label, label_fmt=label_fmt)

    def show_mpl(self, im, enhance=True, clear_fig=True):
        if self._pylab is None:
            import pylab
            self._pylab = pylab
        if self._render_figure is None:
            self._render_figure = self._pylab.figure(1)
        if clear_fig: self._render_figure.clf()

        if enhance:
            nz = im[im > 0.0]
            nim = im / (nz.mean() + 6.0 * np.std(nz))
            nim[nim > 1.0] = 1.0
            nim[nim < 0.0] = 0.0
            del nz
        else:
            nim = im
        ax = self._pylab.imshow(nim[:,:,:3]/nim[:,:,:3].max(), origin='upper')
        return ax

    def draw(self):
        self._pylab.draw()
    
    def save_annotated(self, fn, image, enhance=True, dpi=100, clear_fig=True, 
                       label_fmt=None):
        """
        Save an image with the transfer function represented as a colorbar.

        Parameters
        ----------
        fn : str
           The output filename
        image : ImageArray
           The image to annotate
        enhance : bool, optional
           Enhance the contrast (default: True)
        dpi : int, optional
           Dots per inch in the output image (default: 100)
        clear_fig : bool, optional
           Reset the figure (through pylab.clf()) before drawing.  Setting 
           this to false can allow us to overlay the image onto an 
           existing figure
        label_fmt : str, optional
           A format specifier (e.g., label_fmt="%.2g") to use in formatting 
           the data values that label the transfer function colorbar. 
        
        """
        image = image.swapaxes(0,1) 
        ax = self.show_mpl(image, enhance=enhance, clear_fig=clear_fig)
        self.annotate(ax.axes, enhance, label_fmt=label_fmt)
        self._pylab.savefig(fn, bbox_inches='tight', facecolor='black', dpi=dpi)
        
    def save_image(self, image, fn=None, clip_ratio=None, transparent=False):
        if self.comm.rank == 0 and fn is not None:
            if transparent:
                image.write_png(fn, clip_ratio=clip_ratio, rescale=True,
                                background=None)
            else:
                image.write_png(fn, clip_ratio=clip_ratio, rescale=True,
                                background='black')

    def initialize_source(self):
        return self.volume.initialize_source(self.fields, self.log_fields,
                                             self.no_ghost)

    def get_information(self):
        info_dict = {'fields':self.fields,
                     'type':self.__class__.__name__,
                     'east_vector':self.orienter.unit_vectors[0],
                     'north_vector':self.orienter.unit_vectors[1],
                     'normal_vector':self.orienter.unit_vectors[2],
                     'width':self.width,
                     'dataset':self.ds.fullpath}
        return info_dict

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0, transparent=False):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Parameters
        ----------
        fn : string, optional
            If supplied, the image will be saved out to this before being
            returned.  Scaling will be to the maximum value.
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()
        double_check : bool, optional
            Optionally makes sure that the data contains only valid entries.
            Used for debugging.
        num_threads : int, optional
            If supplied, will use 'num_threads' number of OpenMP threads during
            the rendering.  Defaults to 0, which uses the environment variable
            OMP_NUM_THREADS.
        transparent: bool, optional
            Optionally saves out the 4-channel rgba image, which can appear 
            empty if the alpha channel is low everywhere. Default: False

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
        if num_threads is None:
            num_threads=get_num_threads()
        image = self.new_image()
        args, kwargs = self.get_sampler_args(image)
        sampler = self.get_sampler(args, kwargs)
        self.initialize_source()
        image = ImageArray(self._render(double_check, num_threads, 
                                        image, sampler),
                           info=self.get_information())

        # flip it up/down to handle how the png orientation is done
        image = image[:,::-1,:]
        self.save_image(image, fn=fn, clip_ratio=clip_ratio, 
                       transparent=transparent)
        return image

    def show(self, clip_ratio = None):
        r"""This will take a snapshot and display the resultant image in the
        IPython notebook.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will snapshot and send the resultant
        image to the IPython notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Parameters
        ----------
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Examples
        --------

        >>> cam.show()

        """
        if "__IPYTHON__" in dir(builtins):
            from IPython.core.displaypub import publish_display_data
            image = self.snapshot()[:,:,:3]
            if clip_ratio is not None: clip_ratio *= image.std()
            data = write_bitmap(image, None, clip_ratio)
            publish_display_data(
                data={'image/png': data},
                source='yt.visualization.volume_rendering.camera.Camera',
            )
        else:
            raise YTNotInsideNotebook


    def set_default_light_dir(self):
        self.light_dir = [1.,1.,1.]

    def set_default_light_rgba(self):
        self.light_rgba = [1.,1.,1.,1.]

    def zoom(self, factor):
        r"""Change the distance to the focal point.

        This will zoom the camera in by some `factor` toward the focal point,
        along the current view direction, modifying the left/right and up/down
        extents as well.

        Parameters
        ----------
        factor : float
            The factor by which to reduce the distance to the focal point.


        Notes
        -----

        You will need to call snapshot() again to get a new image.

        """
        self.width /= factor
        self._setup_box_properties(self.width, self.center, self.orienter.unit_vectors)

    def zoomin(self, final, n_steps, clip_ratio = None):
        r"""Loop over a zoomin and return snapshots along the way.

        This will yield `n_steps` snapshots until the current view has been
        zooming in to a final factor of `final`.

        Parameters
        ----------
        final : float
            The zoom factor, with respect to current, desired at the end of the
            sequence.
        n_steps : int
            The number of zoom snapshots to make.
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()


        Examples
        --------

        >>> for i, snapshot in enumerate(cam.zoomin(100.0, 10)):
        ...     iw.write_bitmap(snapshot, "zoom_%04i.png" % i)
        """
        f = final**(1.0/n_steps)
        for i in range(n_steps):
            self.zoom(f)
            yield self.snapshot(clip_ratio = clip_ratio)

    def move_to(self, final, n_steps, final_width=None, exponential=False, clip_ratio = None):
        r"""Loop over a look_at

        This will yield `n_steps` snapshots until the current view has been
        moved to a final center of `final` with a final width of final_width.

        Parameters
        ----------
        final : array_like
            The final center to move to after `n_steps`
        n_steps : int
            The number of look_at snapshots to make.
        final_width: float or array_like, optional
            Specifies the final width after `n_steps`.  Useful for
            moving and zooming at the same time.
        exponential : boolean
            Specifies whether the move/zoom transition follows an
            exponential path toward the destination or linear
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()
            
        Examples
        --------

        >>> for i, snapshot in enumerate(cam.move_to([0.2,0.3,0.6], 10)):
        ...     iw.write_bitmap(snapshot, "move_%04i.png" % i)
        """
        dW = None
        if not isinstance(final, YTArray):
            final = self.ds.arr(final, input_units = "code_length")
        if exponential:
            if final_width is not None:
                if not iterable(final_width):
                    final_width = [final_width, final_width, final_width] 
                if not isinstance(final_width, YTArray):
                    final_width = self.ds.arr(final_width, input_units="code_length")
                    # left/right, top/bottom, front/back 
                if (self.center == 0.0).all():
                    self.center += (final - self.center) / (10. * n_steps)
                final_zoom = final_width/self.width
                dW = final_zoom**(1.0/n_steps)
            else:
                dW = self.ds.arr([1.0,1.0,1.0], "code_length")
            position_diff = final/self.center
            dx = position_diff**(1.0/n_steps)
        else:
            if final_width is not None:
                if not iterable(final_width):
                    final_width = [final_width, final_width, final_width] 
                if not isinstance(final_width, YTArray):
                    final_width = self.ds.arr(final_width, input_units="code_length")
                    # left/right, top/bottom, front/back
                dW = (1.0*final_width-self.width)/n_steps
            else:
                dW = self.ds.arr([0.0,0.0,0.0], "code_length")
            dx = (final-self.center)*1.0/n_steps
        for i in range(n_steps):
            if exponential:
                self.switch_view(center=self.center*dx, width=self.width*dW)
            else:
                self.switch_view(center=self.center+dx, width=self.width+dW)
            yield self.snapshot(clip_ratio = clip_ratio)

    def rotate(self, theta, rot_vector=None):
        r"""Rotate by a given angle

        Rotate the view.  If `rot_vector` is None, rotation will occur
        around the `north_vector`.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to rotate the view.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around
            `north_vector`

        Examples
        --------

        >>> cam.rotate(np.pi/4)
        """
        rotate_all = rot_vector is not None
        if rot_vector is None:
            rot_vector = self.rotation_vector
        else:
            rot_vector = ensure_numpy_array(rot_vector)
            rot_vector = rot_vector/np.linalg.norm(rot_vector)
          
        R = get_rotation_matrix(theta, rot_vector)

        normal_vector = self.front_center-self.center
        normal_vector = normal_vector/np.sqrt((normal_vector**2).sum())

        if rotate_all:
            self.switch_view(
                normal_vector=np.dot(R, normal_vector),
                north_vector=np.dot(R, self.orienter.unit_vectors[1]))
        else:
            self.switch_view(normal_vector=np.dot(R, normal_vector))


    def pitch(self, theta):
        r"""Rotate by a given angle about the horizontal axis

        Pitch the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to pitch the view.

        Examples
        --------

        >>> cam.pitch(np.pi/4)
        """
        rot_vector = self.orienter.unit_vectors[0]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
                normal_vector=np.dot(R, self.orienter.unit_vectors[2]),
                north_vector=np.dot(R, self.orienter.unit_vectors[1]))
        if self.orienter.steady_north:
            self.orienter.north_vector = self.orienter.unit_vectors[1]
 
    def yaw(self, theta):
        r"""Rotate by a given angle about the vertical axis

        Yaw the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to yaw the view.

        Examples
        --------

        >>> cam.yaw(np.pi/4)
        """
        rot_vector = self.orienter.unit_vectors[1]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
                normal_vector=np.dot(R, self.orienter.unit_vectors[2]))
 
    def roll(self, theta):
        r"""Rotate by a given angle about the view normal axis

        Roll the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to roll the view.

        Examples
        --------

        >>> cam.roll(np.pi/4)
        """
        rot_vector = self.orienter.unit_vectors[2]
        R = get_rotation_matrix(theta, rot_vector)
        self.switch_view(
                normal_vector=np.dot(R, self.orienter.unit_vectors[2]),
                north_vector=np.dot(R, self.orienter.unit_vectors[1]))
        if self.orienter.steady_north:
            self.orienter.north_vector = np.dot(R, self.orienter.north_vector)

    def rotation(self, theta, n_steps, rot_vector=None, clip_ratio = None):
        r"""Loop over rotate, creating a rotation

        This will yield `n_steps` snapshots until the current view has been
        rotated by an angle `theta`

        Parameters
        ----------
        theta : float, in radians
            Angle (in radians) by which to rotate the view.
        n_steps : int
            The number of look_at snapshots to make.
        rot_vector  : array_like, optional
            Specify the rotation vector around which rotation will
            occur.  Defaults to None, which sets rotation around the
            original `north_vector`
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Examples
        --------

        >>> for i, snapshot in enumerate(cam.rotation(np.pi, 10)):
        ...     iw.write_bitmap(snapshot, 'rotation_%04i.png' % i)
        """

        dtheta = (1.0*theta)/n_steps
        for i in range(n_steps):
            self.rotate(dtheta, rot_vector=rot_vector)
            yield self.snapshot(clip_ratio = clip_ratio)

data_object_registry["camera"] = Camera

class InteractiveCamera(Camera):
    frames = []

    def snapshot(self, fn = None, clip_ratio = None):
        import matplotlib.pylab as pylab
        pylab.figure(2)
        self.transfer_function.show()
        pylab.draw()
        im = Camera.snapshot(self, fn, clip_ratio)
        pylab.figure(1)
        pylab.imshow(im / im.max())
        pylab.draw()
        self.frames.append(im)

    def rotation(self, theta, n_steps, rot_vector=None):
        for frame in Camera.rotation(self, theta, n_steps, rot_vector):
            if frame is not None:
                self.frames.append(frame)

    def zoomin(self, final, n_steps):
        for frame in Camera.zoomin(self, final, n_steps):
            if frame is not None:
                self.frames.append(frame)

    def clear_frames(self):
        del self.frames
        self.frames = []

    def save(self,fn):
        self._pylab.savefig(fn, bbox_inches='tight', facecolor='black')

    def save_frames(self, basename, clip_ratio=None):
        for i, frame in enumerate(self.frames):
            fn = basename + '_%04i.png'%i
            if clip_ratio is not None:
                write_bitmap(frame, fn, clip_ratio*frame.std())
            else:
                write_bitmap(frame, fn)

data_object_registry["interactive_camera"] = InteractiveCamera

class PerspectiveCamera(Camera):
    r"""A viewpoint into a volume, for perspective volume rendering.

    The camera represents the eye of an observer, which will be used to
    generate ray-cast volume renderings of the domain. The rays start from
    the camera and end on the image plane, which generates a perspective
    view.

    Note: at the moment, this results in a left-handed coordinate
    system view

    Parameters
    ----------
    center : array_like
        The location of the camera
    normal_vector : array_like
        The vector from the camera position to the center of the image plane
    width : float or list of floats
        width[0] and width[1] give the width and height of the image plane, and
        width[2] gives the depth of the image plane (distance between the camera
        and the center of the image plane).
        The view angles thus become:
        2 * arctan(0.5 * width[0] / width[2]) in horizontal direction
        2 * arctan(0.5 * width[1] / width[2]) in vertical direction
    (The following parameters are identical with the definitions in Camera class)
    resolution : int or list of ints
        The number of pixels in each direction.
    transfer_function : `yt.visualization.volume_rendering.TransferFunction`
        The transfer function used to map values to colors in an image.  If
        not specified, defaults to a ProjectionTransferFunction.
    north_vector : array_like, optional
        The 'up' direction for the plane of rays.  If not specific, calculated
        automatically.
    steady_north : bool, optional
        Boolean to control whether to normalize the north_vector
        by subtracting off the dot product of it and the normal
        vector.  Makes it easier to do rotations along a single
        axis.  If north_vector is specified, is switched to
        True. Default: False
    volume : `yt.extensions.volume_rendering.AMRKDTree`, optional
        The volume to ray cast through.  Can be specified for finer-grained
        control, but otherwise will be automatically generated.
    fields : list of fields, optional
        This is the list of fields we want to volume render; defaults to
        Density.
    log_fields : list of bool, optional
        Whether we should take the log of the fields before supplying them to
        the volume rendering mechanism.
    sub_samples : int, optional
        The number of samples to take inside every cell per ray.
    ds : `~yt.data_objects.api.Dataset`
        For now, this is a require parameter!  But in the future it will become
        optional.  This is the dataset to volume render.
    use_kd: bool, optional
        Specifies whether or not to use a kd-Tree framework for
        the Homogenized Volume and ray-casting.  Default to True.
    max_level: int, optional
        Specifies the maximum level to be rendered.  Also
        specifies the maximum level used in the kd-Tree
        construction.  Defaults to None (all levels), and only
        applies if use_kd=True.
    no_ghost: bool, optional
        Optimization option.  If True, homogenized bricks will
        extrapolate out from grid instead of interpolating from
        ghost zones that have to first be calculated.  This can
        lead to large speed improvements, but at a loss of
        accuracy/smoothness in resulting image.  The effects are
        less notable when the transfer function is smooth and
        broad. Default: True
    data_source: data container, optional
        Optionally specify an arbitrary data source to the volume rendering.
        All cells not included in the data source will be ignored during ray
        casting. By default this will get set to ds.all_data().

    """
    def __init__(self, *args, **kwargs):
        Camera.__init__(self, *args, **kwargs)

    def get_sampler_args(self, image):
        east_vec = self.orienter.unit_vectors[0].reshape(3,1)
        north_vec = self.orienter.unit_vectors[1].reshape(3,1)

        px = np.mat(np.linspace(-.5, .5, self.resolution[0]))
        py = np.mat(np.linspace(-.5, .5, self.resolution[1]))

        sample_x = self.width[0] * np.array(east_vec * px).transpose()
        sample_y = self.width[1] * np.array(north_vec * py).transpose()

        vectors = np.zeros((self.resolution[0], self.resolution[1], 3),
                           dtype='float64', order='C')
    
        sample_x = np.repeat(sample_x.reshape(self.resolution[0],1,3), \
                             self.resolution[1], axis=1)
        sample_y = np.repeat(sample_y.reshape(1,self.resolution[1],3), \
                             self.resolution[0], axis=0)

        normal_vec = np.zeros((self.resolution[0], self.resolution[1], 3),
                              dtype='float64', order='C')
        normal_vec[:,:,0] = self.orienter.unit_vectors[2,0]
        normal_vec[:,:,1] = self.orienter.unit_vectors[2,1]
        normal_vec[:,:,2] = self.orienter.unit_vectors[2,2]

        vectors = sample_x + sample_y + normal_vec * self.width[2]

        positions = np.zeros((self.resolution[0], self.resolution[1], 3),
                             dtype='float64', order='C')
        positions[:,:,0] = self.center[0]
        positions[:,:,1] = self.center[1]
        positions[:,:,2] = self.center[2]

        positions = self.ds.arr(positions, input_units="code_length")

        dummy = np.ones(3, dtype='float64')
        image.shape = (self.resolution[0], self.resolution[1],4)

        args = (positions, vectors, self.back_center,
                (0.0,1.0,0.0,1.0),
                image, dummy, dummy,
                np.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args, {'lens_type': 'perspective'}

    def _render(self, double_check, num_threads, image, sampler):
        ncells = sum(b.source_mask.size for b in self.volume.bricks)
        pbar = get_pbar("Ray casting", ncells)
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(self.front_center):
            sampler(brick, num_threads=num_threads)
            total_cells += brick.source_mask.size
            pbar.update(total_cells)

        pbar.finish()
        image = self.finalize_image(sampler.aimage)
        return image

    def finalize_image(self, image):
        view_pos = self.front_center
        image.shape = self.resolution[0], self.resolution[1], 4
        image = self.volume.reduce_tree_images(image, view_pos)
        if self.transfer_function.grey_opacity is False:
            image[:,:,3]=1.0
        return image

    def project_to_plane(self, pos, res=None):
        if res is None:
            res = self.resolution
        sight_vector = pos - self.center
        pos1 = sight_vector
        for i in range(0, sight_vector.shape[0]):
            sight_vector_norm = np.sqrt(np.dot(sight_vector[i], sight_vector[i]))
            sight_vector[i] = sight_vector[i]/sight_vector_norm
        sight_vector = self.ds.arr(sight_vector.value, input_units='dimensionless')
        sight_center = self.center + self.width[2] * self.orienter.unit_vectors[2]

        for i in range(0, sight_vector.shape[0]):
            sight_angle_cos = np.dot(sight_vector[i], self.orienter.unit_vectors[2])
            if np.arccos(sight_angle_cos) < 0.5 * np.pi:
                sight_length = self.width[2] / sight_angle_cos
            else:
                # The corner is on the backwards, then put it outside of the
                # image It can not be simply removed because it may connect to
                # other corner within the image, which produces visible domian
                # boundary line
                sight_length = np.sqrt(self.width[0]**2+self.width[1]**2) / \
                               np.sqrt(1 - sight_angle_cos**2)
            pos1[i] = self.center + sight_length * sight_vector[i]

        dx = np.dot(pos1 - sight_center, self.orienter.unit_vectors[0])
        dy = np.dot(pos1 - sight_center, self.orienter.unit_vectors[1])
        dz = np.dot(pos1 - sight_center, self.orienter.unit_vectors[2])
        # Transpose into image coords.
        px = (res[0]*0.5 + res[0]/self.width[0]*dx).astype('int')
        py = (res[1]*0.5 + res[1]/self.width[1]*dy).astype('int')
        return px, py, dz

    def yaw(self, theta, rot_center):
        r"""Rotate by a given angle about the vertical axis through the
        point center.  This is accomplished by rotating the 
        focal point and then setting the looking vector to point
        to the center.

        Yaw the view.

        Parameters
        ----------
        theta : float, in radians
             Angle (in radians) by which to yaw the view.

        center : a tuple (x, y, z) 
             The point to rotate about

        Examples
        --------

        >>> cam.yaw(np.pi/4, (0., 0., 0.))
        """

        rot_vector = self.orienter.unit_vectors[1]

        focal_point = self.center - rot_center
        R = get_rotation_matrix(theta, rot_vector)
        focal_point = np.dot(R, focal_point) + rot_center

        normal_vector = rot_center - focal_point
        normal_vector = normal_vector/np.sqrt((normal_vector**2).sum())

        self.switch_view(normal_vector=normal_vector, center=focal_point)

    
data_object_registry["perspective_camera"] = PerspectiveCamera

def corners(left_edge, right_edge):
    return np.array([
      [left_edge[:,0], left_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], left_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], right_edge[:,1], left_edge[:,2]],
      [right_edge[:,0], right_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], right_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], left_edge[:,1], right_edge[:,2]],
      [right_edge[:,0], left_edge[:,1], right_edge[:,2]],
      [left_edge[:,0], right_edge[:,1], left_edge[:,2]],
    ], dtype='float64')


class HEALpixCamera(Camera):

    _sampler_object = None 
    
    def __init__(self, center, radius, nside,
                 transfer_function = None, fields = None,
                 sub_samples = 5, log_fields = None, volume = None,
                 ds = None, use_kd=True, no_ghost=False, use_light=False,
                 inner_radius = 10):
        mylog.error('I am sorry, HEALpix Camera does not work yet in 3.0')
        raise NotImplementedError
        ParallelAnalysisInterface.__init__(self)
        if ds is not None: self.ds = ds
        self.center = np.array(center, dtype='float64')
        self.radius = radius
        self.inner_radius = inner_radius
        self.nside = nside
        self.use_kd = use_kd
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function

        if isinstance(self.transfer_function, ProjectionTransferFunction):
            self._sampler_object = InterpolatedProjectionSampler
            self._needs_tf = 0
        else:
            self._sampler_object = VolumeRenderSampler
            self._needs_tf = 1

        if fields is None: fields = ["density"]
        self.fields = fields
        self.sub_samples = sub_samples
        self.log_fields = log_fields
        dd = ds.all_data()
        efields = dd._determine_fields(self.fields)
        if self.log_fields is None:
            self.log_fields = [self.ds._get_field_info(*f).take_log for f in efields]
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if volume is None:
            volume = AMRKDTree(self.ds, data_source=self.data_source)
        self.use_kd = isinstance(volume, AMRKDTree)
        self.volume = volume

    def new_image(self):
        image = np.zeros((12 * self.nside ** 2, 1, 4), dtype='float64', order='C')
        return image

    def get_sampler_args(self, image):
        nv = 12 * self.nside ** 2
        vs = arr_pix2vec_nest(self.nside, np.arange(nv))
        vs.shape = (nv, 1, 3)
        vs += 1e-8
        uv = np.ones(3, dtype='float64')
        positions = np.ones((nv, 1, 3), dtype='float64') * self.center
        dx = min(g.dds.min() for g in self.ds.index.find_point(self.center)[0])
        positions += self.inner_radius * dx * vs
        vs *= self.radius
        args = (positions, vs, self.center,
                (0.0, 1.0, 0.0, 1.0),
                image, uv, uv,
                np.zeros(3, dtype='float64'))
        if self._needs_tf:
            args += (self.transfer_function,)
        args += (self.sub_samples,)
        return args, {}

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting", (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError
        
        view_pos = self.center
        for brick in self.volume.traverse(view_pos):
            sampler(brick, num_threads=num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        
        pbar.finish()
        image = sampler.aimage

        self.finalize_image(image)

        return image

    def finalize_image(self, image):
        view_pos = self.center
        image = self.volume.reduce_tree_images(image, view_pos)
        return image

    def get_information(self):
        info_dict = {'fields':self.fields,
                     'type':self.__class__.__name__,
                     'center':self.center,
                     'radius':self.radius,
                     'dataset':self.ds.fullpath}
        return info_dict


    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0, clim = None, label = None):
        r"""Ray-cast the camera.

        This method instructs the camera to take a snapshot -- i.e., call the ray
        caster -- based on its current settings.

        Parameters
        ----------
        fn : string, optional
            If supplied, the image will be saved out to this before being
            returned.  Scaling will be to the maximum value.
        clip_ratio : float, optional
            If supplied, the 'max_val' argument to write_bitmap will be handed
            clip_ratio * image.std()

        Returns
        -------
        image : array
            An (N,M,3) array of the final returned values, in float64 form.
        """
        if num_threads is None:
            num_threads=get_num_threads()
        image = self.new_image()
        args, kwargs = self.get_sampler_args(image)
        sampler = self.get_sampler(args, kwargs)
        self.volume.initialize_source()
        image = ImageArray(self._render(double_check, num_threads, 
                                        image, sampler),
                           info=self.get_information())
        self.save_image(image, fn=fn, clim=clim, label = label)
        return image

    def save_image(self, image, fn=None, clim=None, label = None):
        if self.comm.rank == 0 and fn is not None:
            # This assumes Density; this is a relatively safe assumption.
            if label is None:
                label = "Projected %s" % (self.fields[0])
            if clim is not None:
                cmin, cmax = clim
            else:
                cmin = cmax = None
            plot_allsky_healpix(image[:,0,0], self.nside, fn, label, 
                                cmin = cmin, cmax = cmax)


class StereoPairCamera(Camera):
    def __init__(self, original_camera, relative_separation = 0.005):
        ParallelAnalysisInterface.__init__(self)
        self.original_camera = original_camera
        self.relative_separation = relative_separation

    def split(self):
        oc = self.original_camera
        uv = oc.orienter.unit_vectors
        c = oc.center
        fc = oc.front_center
        wx, wy, wz = oc.width
        left_normal = fc + uv[1] * 0.5*self.relative_separation * wx - c
        right_normal = fc - uv[1] * 0.5*self.relative_separation * wx - c
        left_camera = Camera(c, left_normal, oc.width,
                             oc.resolution, oc.transfer_function, north_vector=uv[0],
                             volume=oc.volume, fields=oc.fields, log_fields=oc.log_fields,
                             sub_samples=oc.sub_samples, ds=oc.ds)
        right_camera = Camera(c, right_normal, oc.width,
                             oc.resolution, oc.transfer_function, north_vector=uv[0],
                             volume=oc.volume, fields=oc.fields, log_fields=oc.log_fields,
                             sub_samples=oc.sub_samples, ds=oc.ds)
        return (left_camera, right_camera)

class FisheyeCamera(Camera):
    def __init__(self, center, radius, fov, resolution,
                 transfer_function = None, fields = None,
                 sub_samples = 5, log_fields = None, volume = None,
                 ds = None, no_ghost=False, rotation = None, use_light=False):
        ParallelAnalysisInterface.__init__(self)
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        if rotation is None: rotation = np.eye(3)
        self.rotation_matrix = rotation
        self.no_ghost = no_ghost
        if ds is not None: self.ds = ds
        self.center = np.array(center, dtype='float64')
        self.radius = radius
        self.fov = fov
        if iterable(resolution):
            raise RuntimeError("Resolution must be a single int")
        self.resolution = resolution
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        if fields is None: fields = ["density"]
        dd = self.ds.all_data()
        fields = dd._determine_fields(fields)
        self.fields = fields
        if log_fields is None:
            log_fields = [self.ds._get_field_info(*f).take_log for f in fields]
        self.log_fields = log_fields
        self.sub_samples = sub_samples
        if volume is None:
            volume = AMRKDTree(self.ds)
            volume.set_fields(fields, log_fields, no_ghost)
        self.volume = volume

    def get_information(self):
        return {}

    def new_image(self):
        image = np.zeros((self.resolution**2,1,4), dtype='float64', order='C')
        return image
        
    def get_sampler_args(self, image):
        vp = arr_fisheye_vectors(self.resolution, self.fov)
        vp.shape = (self.resolution**2,1,3)
        vp2 = vp.copy()
        for i in range(3):
            vp[:,:,i] = (vp2 * self.rotation_matrix[:,i]).sum(axis=2)
        del vp2
        vp *= self.radius
        uv = np.ones(3, dtype='float64')
        positions = np.ones((self.resolution**2, 1, 3), dtype='float64') * self.center

        args = (positions, vp, self.center,
                (0.0, 1.0, 0.0, 1.0),
                image, uv, uv,
                np.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args, {}


    def finalize_image(self, image):
        image.shape = self.resolution, self.resolution, 4

    def _render(self, double_check, num_threads, image, sampler):
        pbar = get_pbar("Ray casting", (self.volume.brick_dimensions + 1).prod(axis=-1).sum())
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError
        
        view_pos = self.center
        for brick in self.volume.traverse(view_pos):
            sampler(brick, num_threads=num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
            pbar.update(total_cells)
        
        pbar.finish()
        image = sampler.aimage

        self.finalize_image(image)

        return image

class MosaicCamera(Camera):
    def __init__(self, center, normal_vector, width,
                 resolution, transfer_function = None,
                 north_vector = None, steady_north=False,
                 volume = None, fields = None,
                 log_fields = None,
                 sub_samples = 5, ds = None,
                 use_kd=True, l_max=None, no_ghost=True,
                 tree_type='domain',expand_factor=1.0,
                 le=None, re=None, nimx=1, nimy=1, procs_per_wg=None,
                 preload=True, use_light=False):

        ParallelAnalysisInterface.__init__(self)

        self.procs_per_wg = procs_per_wg
        if ds is not None: self.ds = ds
        if not iterable(resolution):
            resolution = (int(resolution/nimx), int(resolution/nimy))
        self.resolution = resolution
        self.nimx = nimx
        self.nimy = nimy
        self.sub_samples = sub_samples
        if not iterable(width):
            width = (width, width, width) # front/back, left/right, top/bottom
        self.width = np.array([width[0], width[1], width[2]])
        self.center = center
        self.steady_north = steady_north
        self.expand_factor = expand_factor
        # This seems to be necessary for now.  Not sure what goes wrong when not true.
        if north_vector is not None: self.steady_north=True
        self.north_vector = north_vector
        self.normal_vector = normal_vector
        if fields is None: fields = ["density"]
        self.fields = fields
        if transfer_function is None:
            transfer_function = ProjectionTransferFunction()
        self.transfer_function = transfer_function
        self.log_fields = log_fields
        self.use_kd = use_kd
        self.l_max = l_max
        self.no_ghost = no_ghost
        self.preload = preload
        
        self.use_light = use_light
        self.light_dir = None
        self.light_rgba = None
        self.le = le
        self.re = re
        self.width[0]/=self.nimx
        self.width[1]/=self.nimy
        
        self.orienter = Orientation(normal_vector, north_vector=north_vector, steady_north=steady_north)
        self.rotation_vector = self.orienter.north_vector
        # self._setup_box_properties(width, center, self.orienter.unit_vectors)
        
        if self.no_ghost:
            mylog.info('Warning: no_ghost is currently True (default). This may lead to artifacts at grid boundaries.')
        self.tree_type = tree_type
        self.volume = volume

        # self.cameras = np.empty(self.nimx*self.nimy)

    def build_volume(self, volume, fields, log_fields, l_max, no_ghost, tree_type, le, re):
        if volume is None:
            if self.use_kd: raise NotImplementedError
            volume = AMRKDTree(self.ds, l_max=l_max, fields=self.fields, 
                               no_ghost=no_ghost, tree_type=tree_type, 
                               log_fields=log_fields, le=le, re=re)
        else:
            self.use_kd = isinstance(volume, AMRKDTree)
        return volume

    def new_image(self):
        image = np.zeros((self.resolution[0], self.resolution[1], 4), dtype='float64', order='C')
        return image

    def _setup_box_properties(self, width, center, unit_vectors):
        owidth = deepcopy(width)
        self.width = width
        self.origin = self.center - 0.5*self.nimx*self.width[0]*self.orienter.unit_vectors[0] \
                                  - 0.5*self.nimy*self.width[1]*self.orienter.unit_vectors[1] \
                                  - 0.5*self.width[2]*self.orienter.unit_vectors[2]
        dx = self.width[0]
        dy = self.width[1]
        offi = (self.imi + 0.5)
        offj = (self.imj + 0.5)
        mylog.info("Mosaic offset: %f %f" % (offi,offj))
        global_center = self.center
        self.center = self.origin
        self.center += offi*dx*self.orienter.unit_vectors[0]
        self.center += offj*dy*self.orienter.unit_vectors[1]
        
        self.box_vectors = np.array([self.orienter.unit_vectors[0]*dx*self.nimx,
                                     self.orienter.unit_vectors[1]*dy*self.nimy,
                                     self.orienter.unit_vectors[2]*self.width[2]])
        self.back_center = self.center - 0.5*self.width[0]*self.orienter.unit_vectors[2]
        self.front_center = self.center + 0.5*self.width[0]*self.orienter.unit_vectors[2]
        self.center = global_center
        self.width = owidth

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0):

        my_storage = {}
        offx,offy = np.meshgrid(range(self.nimx),range(self.nimy))
        offxy = zip(offx.ravel(), offy.ravel())

        for sto, xy in parallel_objects(offxy, self.procs_per_wg, storage = my_storage, 
                                        dynamic=True):
            self.volume = self.build_volume(self.volume, self.fields, self.log_fields, 
                                   self.l_max, self.no_ghost, 
                                   self.tree_type, self.le, self.re)
            self.initialize_source()

            self.imi, self.imj = xy
            mylog.debug('Working on: %i %i' % (self.imi, self.imj))
            self._setup_box_properties(self.width, self.center, self.orienter.unit_vectors)
            image = self.new_image()
            args, kwargs = self.get_sampler_args(image)
            sampler = self.get_sampler(args, kwargs)
            image = self._render(double_check, num_threads, image, sampler)
            sto.id = self.imj*self.nimx + self.imi
            sto.result = image
        image = self.reduce_images(my_storage)
        self.save_image(image, fn=fn, clip_ratio=clip_ratio)
        return image

    def reduce_images(self,im_dict):
        final_image = 0
        if self.comm.rank == 0:
            offx,offy = np.meshgrid(range(self.nimx),range(self.nimy))
            offxy = zip(offx.ravel(), offy.ravel())
            nx,ny = self.resolution
            final_image = np.empty((nx*self.nimx, ny*self.nimy, 4),
                        dtype='float64',order='C')
            for xy in offxy: 
                i, j = xy
                ind = j*self.nimx+i
                final_image[i*nx:(i+1)*nx, j*ny:(j+1)*ny,:] = im_dict[ind]
        return final_image

data_object_registry["mosaic_camera"] = MosaicCamera

def plot_allsky_healpix(image, nside, fn, label = "", rotation = None,
                        take_log = True, resolution=512, cmin=None, cmax=None):
    import matplotlib.figure
    import matplotlib.backends.backend_agg
    if rotation is None: rotation = np.eye(3).astype("float64")

    img, count = pixelize_healpix(nside, image, resolution, resolution, rotation)

    fig = matplotlib.figure.Figure((10, 5))
    ax = fig.add_subplot(1,1,1,projection='aitoff')
    if take_log: func = np.log10
    else: func = lambda a: a
    implot = ax.imshow(func(img), extent=(-np.pi,np.pi,-np.pi/2,np.pi/2),
                       clip_on=False, aspect=0.5, vmin=cmin, vmax=cmax)
    cb = fig.colorbar(implot, orientation='horizontal')
    cb.set_label(label)
    ax.xaxis.set_ticks(())
    ax.yaxis.set_ticks(())
    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    canvas.print_figure(fn)
    return img, count

class ProjectionCamera(Camera):
    def __init__(self, center, normal_vector, width, resolution,
            field, weight=None, volume=None, no_ghost = False, 
            north_vector=None, ds=None, interpolated=False,
            method="integrate"):

        if not interpolated:
            volume = 1

        self.interpolated = interpolated
        self.field = field
        self.weight = weight
        self.resolution = resolution
        self.method = method

        fields = [field]
        if self.weight is not None:
            # This is a temporary field, which we will remove at the end
            # it is given a unique name to avoid conflicting with other 
            # class instances
            self.weightfield = ("index", "temp_weightfield_%u"%(id(self),))
            def _make_wf(f, w):
                def temp_weightfield(a, b):
                    tr = b[f].astype("float64") * b[w]
                    return b.apply_units(tr, a.units)
                    return tr
                return temp_weightfield
            ds.field_info.add_field(self.weightfield,
                function=_make_wf(self.field, self.weight))
            # Now we have to tell the dataset to add it and to calculate
            # its dependencies..
            deps, _ = ds.field_info.check_derived_fields([self.weightfield])
            ds.field_dependencies.update(deps)
            fields = [self.weightfield, self.weight]
        
        self.fields = fields
        self.log_fields = [False]*len(self.fields)
        Camera.__init__(self, center, normal_vector, width, resolution, None,
                fields = fields, ds=ds, volume=volume,
                log_fields=self.log_fields, 
                north_vector=north_vector,
                no_ghost=no_ghost)

    # this would be better in an __exit__ function, but that would require
    # changes in code that uses this class
    def __del__(self):
        if hasattr(self,"weightfield") and hasattr(self,"ds"):
            try:
                self.ds.field_info.pop(self.weightfield)
                self.ds.field_dependencies.pop(self.weightfield)
            except KeyError:
                pass
        try:
            Camera.__del__(self)
        except AttributeError:
            pass

    def get_sampler(self, args, kwargs):
        if self.interpolated:
            sampler = InterpolatedProjectionSampler(*args, **kwargs)
        else:
            sampler = ProjectionSampler(*args, **kwargs)
        return sampler

    def initialize_source(self):
        if self.interpolated:
            Camera.initialize_source(self)
        else:
            pass

    def get_sampler_args(self, image):
        rotp = np.concatenate([self.orienter.inv_mat.ravel('F'), self.back_center.ravel()])
        args = (np.atleast_3d(rotp), np.atleast_3d(self.box_vectors[2]),
                self.back_center,
            (-self.width[0]/2., self.width[0]/2.,
             -self.width[1]/2., self.width[1]/2.),
            image, self.orienter.unit_vectors[0], self.orienter.unit_vectors[1],
                np.array(self.width, dtype='float64'), self.sub_samples)
        return args, {'lens_type': 'plane-parallel'}

    def finalize_image(self,image):
        ds = self.ds
        dd = ds.all_data()
        field = dd._determine_fields([self.field])[0]
        finfo = ds._get_field_info(*field)
        dl = 1.0
        if self.method == "integrate":
            if self.weight is None:
                dl = self.width[2].in_units(ds.unit_system["length"])
            else:
                image[:, : ,0] /= image[:, :, 1]

        return ImageArray(image[:, :, 0], finfo.units, 
                          registry=ds.unit_registry) * dl


    def _render(self, double_check, num_threads, image, sampler):
        # Calculate the eight corners of the box
        # Back corners ...
        if self.interpolated:
            return Camera._render(self, double_check, num_threads, image,
                    sampler)
        ds = self.ds
        width = self.width[2]
        north_vector = self.orienter.unit_vectors[0]
        east_vector = self.orienter.unit_vectors[1]
        normal_vector = self.orienter.unit_vectors[2]
        fields = self.fields

        mi = ds.domain_right_edge.copy()
        ma = ds.domain_left_edge.copy()
        for off1 in [-1, 1]:
            for off2 in [-1, 1]:
                for off3 in [-1, 1]:
                    this_point = (self.center + width/2. * off1 * north_vector
                                         + width/2. * off2 * east_vector
                                         + width/2. * off3 * normal_vector)
                    np.minimum(mi, this_point, mi)
                    np.maximum(ma, this_point, ma)
        # Now we have a bounding box.
        data_source = ds.region(self.center, mi, ma)

        for i, (grid, mask) in enumerate(data_source.blocks):
            data = [(grid[field] * mask).astype("float64") for field in fields]
            pg = PartitionedGrid(
                grid.id, data,
                mask.astype('uint8'),
                grid.LeftEdge, grid.RightEdge, grid.ActiveDimensions.astype("int64"))
            grid.clear_data()
            sampler(pg, num_threads = num_threads)

        image = self.finalize_image(sampler.aimage)
        return image

    def save_image(self, image, fn=None, clip_ratio=None):
        dd = self.ds.all_data()
        field = dd._determine_fields([self.field])[0]
        finfo = self.ds._get_field_info(*field)
        if finfo.take_log:
            im = np.log10(image)
        else:
            im = image
        if self.comm.rank == 0 and fn is not None:
            if clip_ratio is not None:
                write_image(im, fn)
            else:
                write_image(im, fn)

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0):

        if num_threads is None:
            num_threads=get_num_threads()

        image = self.new_image()

        args, kwargs = self.get_sampler_args(image)

        sampler = self.get_sampler(args, kwargs)

        self.initialize_source()

        image = ImageArray(self._render(double_check, num_threads, 
                                        image, sampler),
                           info=self.get_information())

        self.save_image(image, fn=fn, clip_ratio=clip_ratio)

        return image
    snapshot.__doc__ = Camera.snapshot.__doc__

data_object_registry["projection_camera"] = ProjectionCamera

class SphericalCamera(Camera):
    def __init__(self, *args, **kwargs):
        Camera.__init__(self, *args, **kwargs)
        if(self.resolution[0]/self.resolution[1] != 2):
            mylog.info('Warning: It\'s recommended to set the aspect ratio to 2:1')
        self.resolution = np.asarray(self.resolution) + 2

    def get_sampler_args(self, image):
        px = np.linspace(-np.pi, np.pi, self.resolution[0], endpoint=True)[:,None]
        py = np.linspace(-np.pi/2., np.pi/2., self.resolution[1], endpoint=True)[None,:]
        
        vectors = np.zeros((self.resolution[0], self.resolution[1], 3),
                           dtype='float64', order='C')
        vectors[:,:,0] = np.cos(px) * np.cos(py)
        vectors[:,:,1] = np.sin(px) * np.cos(py)
        vectors[:,:,2] = np.sin(py)

        vectors = vectors * self.width[0]
        positions = self.center + vectors * 0
        R1 = get_rotation_matrix(0.5*np.pi, [1,0,0])
        R2 = get_rotation_matrix(0.5*np.pi, [0,0,1])
        uv = np.dot(R1, self.orienter.unit_vectors)
        uv = np.dot(R2, uv)
        vectors.reshape((self.resolution[0]*self.resolution[1], 3))
        vectors = np.dot(vectors, uv)
        vectors.reshape((self.resolution[0], self.resolution[1], 3))

        dummy = np.ones(3, dtype='float64')
        image.shape = (self.resolution[0]*self.resolution[1],1,4)
        vectors.shape = (self.resolution[0]*self.resolution[1],1,3)
        positions.shape = (self.resolution[0]*self.resolution[1],1,3)
        args = (positions, vectors, self.back_center,
                (0.0,1.0,0.0,1.0),
                image, dummy, dummy,
                np.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args, {'lens_type': 'spherical'}

    def _render(self, double_check, num_threads, image, sampler):
        ncells = sum(b.source_mask.size for b in self.volume.bricks)
        pbar = get_pbar("Ray casting", ncells)
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(self.front_center):
            sampler(brick, num_threads=num_threads)
            total_cells += brick.source_mask.size
            pbar.update(total_cells)

        pbar.finish()
        image = self.finalize_image(sampler.aimage)
        return image

    def finalize_image(self, image):
        view_pos = self.front_center
        image.shape = self.resolution[0], self.resolution[1], 4
        image = self.volume.reduce_tree_images(image, view_pos)
        if self.transfer_function.grey_opacity is False:
            image[:,:,3]=1.0
        image = image[1:-1,1:-1,:]
        return image

data_object_registry["spherical_camera"] = SphericalCamera

class StereoSphericalCamera(Camera):
    def __init__(self, *args, **kwargs):
        self.disparity = kwargs.pop('disparity', 0.)
        Camera.__init__(self, *args, **kwargs)
        self.disparity = self.ds.arr(self.disparity, input_units="code_length")
        self.disparity_s = self.ds.arr(0., input_units="code_length")
        if(self.resolution[0]/self.resolution[1] != 2):
            mylog.info('Warning: It\'s recommended to set the aspect ratio to be 2:1')
        self.resolution = np.asarray(self.resolution) + 2
        if(self.disparity<=0.):
            self.disparity = self.width[0]/1000.
            mylog.info('Warning: Invalid value of disparity; ' \
                       'now reset it to %f' % self.disparity)

    def get_sampler_args(self, image):
        px = np.linspace(-np.pi, np.pi, self.resolution[0], endpoint=True)[:,None]
        py = np.linspace(-np.pi/2., np.pi/2., self.resolution[1], endpoint=True)[None,:]

        vectors = np.zeros((self.resolution[0], self.resolution[1], 3),
                           dtype='float64', order='C')
        vectors[:,:,0] = np.cos(px) * np.cos(py)
        vectors[:,:,1] = np.sin(px) * np.cos(py)
        vectors[:,:,2] = np.sin(py)
        vectors2 = np.zeros((self.resolution[0], self.resolution[1], 3), 
                            dtype='float64', order='C')
        vectors2[:,:,0] = -np.sin(px) * np.ones((1, self.resolution[1]))
        vectors2[:,:,1] = np.cos(px) * np.ones((1, self.resolution[1]))
        vectors2[:,:,2] = 0

        positions = self.center + vectors2 * self.disparity_s
        vectors = vectors * self.width[0]
        R1 = get_rotation_matrix(0.5*np.pi, [1,0,0])
        R2 = get_rotation_matrix(0.5*np.pi, [0,0,1])
        uv = np.dot(R1, self.orienter.unit_vectors)
        uv = np.dot(R2, uv)
        vectors.reshape((self.resolution[0]*self.resolution[1], 3))
        vectors = np.dot(vectors, uv)
        vectors.reshape((self.resolution[0], self.resolution[1], 3))

        dummy = np.ones(3, dtype='float64')
        image.shape = (self.resolution[0]*self.resolution[1],1,4)
        vectors.shape = (self.resolution[0]*self.resolution[1],1,3)
        positions.shape = (self.resolution[0]*self.resolution[1],1,3)
        args = (positions, vectors, self.back_center,
                (0.0,1.0,0.0,1.0),
                image, dummy, dummy,
                np.zeros(3, dtype='float64'),
                self.transfer_function, self.sub_samples)
        return args, {'lens_type': 'stereo-spherical'}

    def snapshot(self, fn = None, clip_ratio = None, double_check = False,
                 num_threads = 0, transparent=False):
        
        if num_threads is None:
            num_threads=get_num_threads()

        self.disparity_s = self.disparity
        image1 = self.new_image()
        args1, kwargs1 = self.get_sampler_args(image1)
        sampler1 = self.get_sampler(args1, kwargs1)
        self.initialize_source()
        image1 = self._render(double_check, num_threads,
                              image1, sampler1, '(Left) ')

        self.disparity_s = -self.disparity
        image2 = self.new_image()
        args2, kwargs2 = self.get_sampler_args(image2)
        sampler2 = self.get_sampler(args2, kwargs2)
        self.initialize_source()
        image2 = self._render(double_check, num_threads,
                              image2, sampler2, '(Right)')

        image = np.hstack([image1, image2])
        image = self.volume.reduce_tree_images(image, self.center)
        image = ImageArray(image, info = self.get_information())
        self.save_image(image, fn=fn, clip_ratio=clip_ratio,
                        transparent=transparent)
        return image

    def _render(self, double_check, num_threads, image, sampler, msg):
        ncells = sum(b.source_mask.size for b in self.volume.bricks)
        pbar = get_pbar("Ray casting "+msg, ncells)
        total_cells = 0
        if double_check:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(self.front_center):
            sampler(brick, num_threads=num_threads)
            total_cells += brick.source_mask.size 
            pbar.update(total_cells)

        pbar.finish()

        image = sampler.aimage.copy()
        image.shape = self.resolution[0], self.resolution[1], 4
        if self.transfer_function.grey_opacity is False:
            image[:,:,3]=1.0
        image = image[1:-1,1:-1,:]
        return image

data_object_registry["stereospherical_camera"] = StereoSphericalCamera

def off_axis_projection(ds, center, normal_vector, width, resolution,
                        field, weight = None, 
                        volume = None, no_ghost = False, interpolated = False,
                        north_vector = None, method = "integrate"):
    r"""Project through a dataset, off-axis, and return the image plane.

    This function will accept the necessary items to integrate through a volume
    at an arbitrary angle and return the integrated field of view to the user.
    Note that if a weight is supplied, it will multiply the pre-interpolated
    values together, then create cell-centered values, then interpolate within
    the cell to conduct the integration.

    Parameters
    ----------
    ds : `~yt.data_objects.api.Dataset`
        This is the dataset to volume render.
    center : array_like
        The current 'center' of the view port -- the focal point for the
        camera.
    normal_vector : array_like
        The vector between the camera position and the center.
    width : float or list of floats
        The current width of the image.  If a single float, the volume is
        cubical, but if not, it is left/right, top/bottom, front/back
    resolution : int or list of ints
        The number of pixels in each direction.
    field : string
        The field to project through the volume
    weight : optional, default None
        If supplied, the field will be pre-multiplied by this, then divided by
        the integrated value of this field.  This returns an average rather
        than a sum.
    volume : `yt.extensions.volume_rendering.AMRKDTree`, optional
        The volume to ray cast through.  Can be specified for finer-grained
        control, but otherwise will be automatically generated.
    no_ghost: bool, optional
        Optimization option.  If True, homogenized bricks will
        extrapolate out from grid instead of interpolating from
        ghost zones that have to first be calculated.  This can
        lead to large speed improvements, but at a loss of
        accuracy/smoothness in resulting image.  The effects are
        less notable when the transfer function is smooth and
        broad. Default: True
    interpolated : optional, default False
        If True, the data is first interpolated to vertex-centered data, 
        then tri-linearly interpolated along the ray. Not suggested for 
        quantitative studies.
    method : string
         The method of projection.  Valid methods are:

         "integrate" with no weight_field specified : integrate the requested
         field along the line of sight.

         "integrate" with a weight_field specified : weight the requested
         field by the weighting field and integrate along the line of sight.

         "sum" : This method is the same as integrate, except that it does not
         multiply by a path length when performing the integration, and is
         just a straight summation of the field along the given axis. WARNING:
         This should only be used for uniform resolution grid datasets, as other
         datasets may result in unphysical images.

    Returns
    -------
    image : array
        An (N,N) array of the final integrated values, in float64 form.

    Examples
    --------

    >>> image = off_axis_projection(ds, [0.5, 0.5, 0.5], [0.2,0.3,0.4],
                      0.2, N, "temperature", "density")
    >>> write_image(np.log10(image), "offaxis.png")

    """
    projcam = ProjectionCamera(center, normal_vector, width, resolution,
                               field, weight=weight, ds=ds, volume=volume,
                               no_ghost=no_ghost, interpolated=interpolated, 
                               north_vector=north_vector, method=method)
    image = projcam.snapshot()
    return image[:,:]

