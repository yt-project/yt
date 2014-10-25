"""
RenderSource Class



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog
from yt.units import yt_array
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from transfer_function_helper import TransferFunctionHelper
from transfer_functions import TransferFunction
from utils import new_volume_render_sampler, data_source_or_all, \
    get_corners

from zbuffer_array import ZBuffer
from yt.utilities.lib.misc_utilities import \
    lines


class RenderSource(ParallelAnalysisInterface):

    """Base Class for Render Sources. Will be inherited for volumes,
       streamlines, etc"""

    def __init__(self):
        mylog.debug("Entering %s" % str(self))
        super(RenderSource, self).__init__()
        self.opaque = False
        self.zbuffer = None

    def render(self, camera, zbuffer=None):
        pass

    def validate(self):
        pass


class OpaqueSource(RenderSource):
    """docstring for OpaqueSource"""
    def __init__(self):
        mylog.debug("Entering %s" % str(self))
        super(OpaqueSource, self).__init__()
        self.opaque = True

    def set_zbuffer(self, zbuffer):
        self.zbuffer = zbuffer

    def render(self, camera, zbuffer=None):
        # This is definitely wrong for now
        return self.zbuffer


class VolumeSource(RenderSource):

    """docstring for VolumeSource"""
    _image = None
    data_source = None

    def __init__(self, data_source, field, auto=True):
        #mylog.debug("Entering %s" % str(self))
        super(VolumeSource, self).__init__()
        self.data_source = data_source_or_all(data_source)
        self.field = field
        self.volume = None
        self.current_image = None
        self.double_check = False
        self.num_threads = 0
        self.num_samples = 10

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
        if not isinstance(transfer_function, TransferFunction):
            raise RuntimeError("transfer_function not of correct type")
        self.transfer_function = transfer_function
        return self

    def validate(self):
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
        log_field = self.data_source.pf.field_info[field].take_log
        self.volume.set_fields([field], [log_field], no_ghost)
        self.field = field

    def set_sampler(self, camera, sampler_type='volume-render'):
        """docstring for add_sampler"""
        if sampler_type == 'volume-render':
            sampler = new_volume_render_sampler(camera, self)
        else:
            NotImplementedError("%s not implemented yet" % sampler_type)
        self.sampler = sampler
        assert(self.sampler is not None)

    def render(self, camera, zbuffer=None):
        camera.lens.set_camera(camera)
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
            self.sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
        mylog.debug("Done casting rays")

        self.current_image = self.finalize_image(camera, self.sampler.aimage)
        if zbuffer is None:
            self.zbuffer = ZBuffer(self.current_image,
                                   np.inf*np.ones(self.current_image.shape[:2]))
        else:
            newz = zbuffer.z.copy()
            newz[:] = 0.0
            self.zbuffer = ZBuffer(self.current_image, newz)

        return self.current_image

    def finalize_image(self, camera, image):
        image.shape = camera.resolution[0], camera.resolution[1], 4
        image = self.volume.reduce_tree_images(image,
                                               camera.lens.viewpoint)
        if self.transfer_function.grey_opacity is False:
            image[:, :, 3] = 1.0
        return image

    def __repr__(self):
        disp = "<Volume Source>:%s" % str(self.data_source)
        return disp


class LineSource(OpaqueSource):

    """Add set of opaque lines to a scene."""
    _image = None
    data_source = None

    def __init__(self, positions, colors=None, color_stride=1):
        super(LineSource, self).__init__()
        self.positions = positions
        if colors is None:
            color = np.ones([positions.shape[0], 4])
        self.colors = colors
        self.color_stride = color_stride

    def render(self, camera, zbuffer=None):
        vertices = self.positions
        empty = camera.lens.new_image(camera)
        z = np.empty(empty.shape[:2])
        empty[:] = 0.0
        z[:] = np.inf
        zbuff = ZBuffer(empty, z)

        # DRAW SOME LINES
        camera.lens.setup_box_properties(camera)
        px, py, dz = camera.lens.project_to_plane(camera, vertices)
        lines(empty, px, py, self.colors, self.color_stride)

        zdummy = -np.ones_like(empty)
        lines(zdummy, px, py, np.vstack([dz, dz, dz, dz]).T, 1)
        z[:, :] = zdummy[:, :, 3]
        z[z == -1] = np.inf
        print 'When rendering opaque lines:', z.min(), z.max()
        self.zbuffer = zbuff
        return zbuff

    def __repr__(self):
        disp = "<Line Source>"
        return disp

class BoxSource(LineSource):
    """Add a box to the scene"""
    def __init__(self, left_edge, right_edge, color=None):
        if color is None:
            color = np.array([1.0,1.0,1.0,1.0]) 
            color.shape = (1, 4)
        corners = get_corners(left_edge.copy(), right_edge.copy())
        order  = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]
        vertices = np.empty([24,3])
        for i in xrange(3):
            vertices[:,i] = corners[order,i,...].ravel(order='F')
        print vertices
        super(BoxSource, self).__init__(vertices, color, color_stride=24)
