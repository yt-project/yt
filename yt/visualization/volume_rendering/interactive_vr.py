# This is for the experimental OpenGL Volume Rendering

import os
import re
import OpenGL.GL as GL
from collections import OrderedDict

import numpy as np
from yt.utilities.math_utils import \
    get_translate_matrix, \
    get_scale_matrix, \
    get_lookat_matrix, \
    get_perspective_matrix, \
    quaternion_mult, \
    quaternion_to_rotation_matrix, \
    rotation_matrix_to_quaternion
from yt.utilities.exceptions import \
    YTInvalidShaderType
from .shader_objects import known_shaders, ShaderProgram

import matplotlib.cm as cm

def _compile_shader(source, shader_type=None):
    if shader_type is None:
        try:
            shader_type = re.match("^.*\.(vertex|fragment)shader$",
                                   source).groups()[0]
        except AttributeError:
            raise YTInvalidShaderType(source)

    if os.path.isfile(source):
        sh_directory = ''
    else:
        sh_directory = os.path.join(os.path.dirname(__file__), "shaders")
    shader = GL.glCreateShader(
        eval("GL.GL_{}_SHADER".format(shader_type.upper()))
    )

    with open(os.path.join(sh_directory, source), 'r') as fp:
        GL.glShaderSource(shader, fp.read())
    GL.glCompileShader(shader)
    result = GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS)
    if not(result):
        raise RuntimeError(GL.glGetShaderInfoLog(shader))
    return shader

def link_shader_program(shaders):
    """Create a shader program with from compiled shaders."""
    compiled_shaders = []
    program = GL.glCreateProgram()
    for shader in shaders:
        compiled_shaders.append(_compile_shader(shader))
        GL.glAttachShader(program, compiled_shaders[-1])
    GL.glLinkProgram(program)
    # check linking error
    result = GL.glGetProgramiv(program, GL.GL_LINK_STATUS)
    if not(result):
        raise RuntimeError(GL.glGetProgramInfoLog(program))
    for shader in compiled_shaders:
        GL.glDeleteShader(shader)
    return program


bbox_vertices = np.array(
      [[ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  1.,  1.],
       [ 0.,  1.,  1.,  1.],
       [ 1.,  1.,  0.,  1.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  1.,  0.,  1.],
       [ 1.,  0.,  1.,  1.],
       [ 0.,  0.,  0.,  1.],
       [ 1.,  0.,  0.,  1.],
       [ 1.,  1.,  0.,  1.],
       [ 1.,  0.,  0.,  1.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  1.,  1.,  1.],
       [ 0.,  1.,  0.,  1.],
       [ 1.,  0.,  1.,  1.],
       [ 0.,  0.,  1.,  1.],
       [ 0.,  0.,  0.,  1.],
       [ 0.,  1.,  1.,  1.],
       [ 0.,  0.,  1.,  1.],
       [ 1.,  0.,  1.,  1.],
       [ 1.,  1.,  1.,  1.],
       [ 1.,  0.,  0.,  1.],
       [ 1.,  1.,  0.,  1.],
       [ 1.,  0.,  0.,  1.],
       [ 1.,  1.,  1.,  1.],
       [ 1.,  0.,  1.,  1.],
       [ 1.,  1.,  1.,  1.],
       [ 1.,  1.,  0.,  1.],
       [ 0.,  1.,  0.,  1.],
       [ 1.,  1.,  1.,  1.],
       [ 0.,  1.,  0.,  1.],
       [ 0.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.],
       [ 0.,  1.,  1.,  1.],
       [ 1.,  0.,  1.,  1.]], dtype=np.float32)

FULLSCREEN_QUAD = np.array(
    [-1.0, -1.0, 0.0,
     +1.0, -1.0, 0.0,
     -1.0, +1.0, 0.0,
     -1.0, +1.0, 0.0,
     +1.0, -1.0, 0.0,
     +1.0, +1.0, 0.0], dtype=np.float32
)

class TrackballCamera(object):
    """

    This class implements a basic "Trackball" or "Arcball" camera control system
    that allows for unconstrained 3D rotations without suffering from Gimbal lock.
    Following Ken Shoemake's orginal C implementation (Graphics Gems IV, III.1)
    we project mouse movements onto the unit sphere and use quaternions to
    represent the corresponding rotation.

    See also:
    https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Arcball

    """

    def __init__(self,
                 position=(0.0, 0.0, 1.0),
                 focus=(0.0, 0.0, 0.0),
                 up=(0.0, 1.0, 0.0),
                 fov=45.0, near_plane=0.01, far_plane=20.0,
                 aspect_ratio=8.0/6.0):
        self.view_matrix = np.zeros((4, 4), dtype=np.float32)
        self.proj_matrix = np.zeros((4, 4), dtype=np.float32)
        self.proj_func = get_perspective_matrix
        self.position = np.array(position)
        self.focus = np.array(focus)
        self.fov = fov
        self.near_plane = near_plane
        self.far_plane = far_plane
        self.aspect_ratio = aspect_ratio
        self.up = np.array(up)
        cmap = cm.get_cmap("algae")
        self.cmap = np.array(cmap(np.linspace(0, 1, 256)), dtype=np.float32)
        self.cmap_min = 1e55
        self.cmap_max = -1e55
        self.cmap_log = True
        self.cmap_new = True

        self.view_matrix = get_lookat_matrix(self.position,
                                             self.focus,
                                             self.up)

        rotation_matrix = self.view_matrix[0:3,0:3]
        self.orientation = rotation_matrix_to_quaternion(rotation_matrix)

    def _map_to_surface(self, mouse_x, mouse_y):
        # right now this just maps to the surface of the unit sphere
        x, y = mouse_x, mouse_y
        mag = np.sqrt(x*x + y*y)
        if (mag > 1.0):
            x /= mag
            y /= mag
            z = 0.0
        else:
            z = np.sqrt(1.0 - mag**2)
        return np.array([x, -y, z])

    def update_orientation(self, start_x, start_y, end_x, end_y):
        old = self._map_to_surface(start_x, start_y)
        new = self._map_to_surface(end_x, end_y)

        # dot product controls the angle of the rotation
        w = old[0]*new[0] + old[1]*new[1] + old[2]*new[2]

        # cross product gives the rotation axis
        x = old[1]*new[2] - old[2]*new[1]
        y = old[2]*new[0] - old[0]*new[2]
        z = old[0]*new[1] - old[1]*new[0]

        q = np.array([w, x, y, z])

        #renormalize to prevent floating point issues
        mag = np.sqrt(w**2 + x**2 + y**2 + z**2)
        q /= mag

        self.orientation = quaternion_mult(self.orientation, q)

    def compute_matrices(self):
        rotation_matrix = quaternion_to_rotation_matrix(self.orientation)
        dp = np.linalg.norm(self.position - self.focus)*rotation_matrix[2]
        self.position = dp + self.focus
        self.up = rotation_matrix[1]

        self.view_matrix = get_lookat_matrix(self.position,
                                             self.focus,
                                             self.up)

        self.projection_matrix = self.proj_func(self.fov,
                                                self.aspect_ratio,
                                                self.near_plane,
                                                self.far_plane)

    def get_viewpoint(self):
        return self.position

    def get_view_matrix(self):
        return self.view_matrix

    def get_projection_matrix(self):
        return self.projection_matrix


    def update_cmap_minmax(self, minval, maxval, iflog):
        self.cmap_log = iflog
        self.cmap_min = minval
        self.cmap_max = maxval

class Camera:
    def __init__(self, position = (0, 0, 0), fov = 60.0, near_plane = 0.01,
            far_plane = 20, aspect_ratio = 8.0 / 6.0, focus = (0, 0, 0),
            up = (0, 0, 1)):
        self.position = np.array(position)
        self.fov = fov
        self.near_plane = near_plane
        self.far_plane = far_plane
        self.aspect_ratio = aspect_ratio
        self.up = np.array(up)
        self.focus = np.array(focus)

    def get_viewpoint(self):
        return self.position

    def get_view_matrix(self):
        return get_lookat_matrix(self.position, self.focus, self.up)

    def get_projection_matrix(self):
        return self.proj_func(self.fov, self.aspect_ratio, self.near_plane,
                              self.far_plane)

    def update_position(self, theta, phi):
        rho = np.linalg.norm(self.position)
        curr_theta = np.arctan2( self.position[1], self.position[0] )
        curr_phi = np.arctan2( np.linalg.norm(self.position[:2]), self.position[2])

        curr_theta += theta
        curr_phi += phi

        self.position[0] = rho * np.sin(curr_phi) * np.cos(curr_theta)
        self.position[1] = rho * np.sin(curr_phi) * np.sin(curr_theta)
        self.position[2] = rho * np.cos(curr_phi)

class SceneComponent(object):
    name = None
    _program = None
    _program_invalid = True
    fragment_shader = None
    vertex_shader = None
    def __init__(self):
        self.vert_attrib = OrderedDict()
        self.vert_arrays = OrderedDict()

    @property
    def program(self):
        if self._program_invalid:
            if self._program is not None:
                self._program.delete_program()
            self._program = ShaderProgram(self.vertex_shader,
                self.fragment_shader)
            self._program_invalid = False
        return self._program

    def _initialize_vertex_array(self, name):
        if name in self.vert_arrays:
            GL.glDeleteVertexArrays(1, [self.vert_arrays[name]])
        self.vert_arrays[name] = GL.glGenVertexArrays(1)

    def run_program(self):
        with self.program.enable():
            if len(self.vert_arrays) != 1:
                raise NotImplementedError
            for vert_name in self.vert_arrays:
                GL.glBindVertexArray(self.vert_arrays[vert_name])
            for an in self.vert_attrib:
                bind_loc, size = self.vert_attrib[an]
                self.program.bind_vert_attrib(an, bind_loc, size)
            self._set_uniforms(self.program)
            self.draw()
            for an in self.vert_attrib:
                self.program.disable_vert_attrib(an)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

    def add_vert_attrib(self, name, arr, each):
        self.vert_attrib[name] = (GL.glGenBuffers(1), each)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vert_attrib[name][0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, arr.nbytes, arr, GL.GL_STATIC_DRAW)

    def set_shader(self, name):
        r""" Compiles and links a fragment shader from a set of known shaders.

        Parameters
        ----------
        name : String
            The name of the fragment shader to use.

        """
        shader = known_shaders[name]()
        if shader.shader_type == "vertex":
            self.vertex_shader = shader
        elif shader.shader_type == "fragment":
            self.fragment_shader = shader
        else:
            raise KeyError(shader.shader_type)
        self._program_invalid = True

class BlockCollection(SceneComponent):
    name = "block_collection"
    def __init__(self):
        super(BlockCollection, self).__init__()
        self.set_shader("default.v")
        self.set_shader("max_intensity.f")
        self.data_source = None

        self.blocks = {} # A collection of PartionedGrid objects
        self.block_order = []

        self.gl_texture_names = []

        self.redraw = True

        self.geometry_loaded = False

        self.camera = None
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthFunc(GL.GL_LESS)

        self._init_blending()

    def _init_blending(self):
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendColor(1.0, 1.0, 1.0, 1.0)
        GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
        GL.glBlendEquation(GL.GL_MAX)

    def set_fields_log(self, log_field):
        self.add_data(self.data_source, self.data_source.tiles.fields[0], log_field)

    def add_data(self, data_source, field, log_field=True):
        r"""Adds a source of data for the block collection.

        Given a `data_source` and a `field` to populate from, adds the data
        to the block collection so that is able to be rendered.

        Parameters
        ----------
        data_source : YTRegion
            A YTRegion object to use as a data source.
        field - String
            A field to populate from.

        """
        self.data_source = data_source
        self.data_source.tiles.set_fields([field], [log_field], no_ghost=False)
        self.data_logged = log_field
        self.blocks = {}
        self.block_order = []
        # Every time we change our data source, we wipe all existing ones.
        # We now set up our vertices into our current data source.
        vert, dx, le, re = [], [], [], []
        self.min_val = 1e60
        self.max_val = -1e60
        for i, block in enumerate(self.data_source.tiles.traverse()):
            self.min_val = min(self.min_val, np.nanmin(block.my_data[0].min()))
            self.max_val = max(self.max_val, np.nanmax(block.my_data[0].max()))
            self.blocks[id(block)] = (i, block)
            vert.append(self._compute_geometry(block, bbox_vertices))
            dds = (block.RightEdge - block.LeftEdge)/block.my_data[0].shape
            n = vert[-1].size/4
            dx.append([dds.astype('f4') for _ in range(n)])
            le.append([block.LeftEdge.astype('f4') for _ in range(n)])
            re.append([block.RightEdge.astype('f4') for _ in range(n)])
            self.block_order.append(id(block))

        LE = np.array([b.LeftEdge
                       for i, b in self.blocks.values()]).min(axis=0)
        RE = np.array([b.RightEdge
                       for i, b in self.blocks.values()]).max(axis=0)
        self.diagonal = np.sqrt(((RE - LE) ** 2).sum())
        # Now we set up our buffer
        vert = np.concatenate(vert)
        dx = np.concatenate(dx)
        le = np.concatenate(le)
        re = np.concatenate(re)

        self._initialize_vertex_array("block_info")
        self.add_vert_attrib("model_vertex", vert, 4)
        self.add_vert_attrib("in_dx", dx, 3)
        self.add_vert_attrib("in_left_edge", le, 3)
        self.add_vert_attrib("in_right_edge", re, 3)

        # Now we set up our textures
        self._load_textures()

    def set_camera(self, camera):
        r"""Sets the camera for the block collection.

        Parameters
        ----------
        camera : Camera
            A simple camera object.

        """
        self.block_order = []
        for block in self.data_source.tiles.traverse(viewpoint = camera.position):
            self.block_order.append(id(block))
        self.camera = camera
        self.redraw = True

    def draw(self):
        r"""Runs a given shader program on the block collection.  It is assumed
        that the GL Context has been set up, which is typically handled by the
        run_program method.
        """

        # clear the color and depth buffer
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glActiveTexture(GL.GL_TEXTURE0)

        for bi in self.block_order:
            tex_i, block = self.blocks[bi]
            ti = self.gl_texture_names[tex_i]
            GL.glBindTexture(GL.GL_TEXTURE_3D, ti)
            GL.glDrawArrays(GL.GL_TRIANGLES, tex_i*36, 36)

    def _set_uniforms(self, shader_program):
        shader_program._set_uniform("projection",
                self.camera.get_projection_matrix())
        shader_program._set_uniform("lookat",
                self.camera.get_view_matrix())
        shader_program._set_uniform("viewport",
                np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype = 'f4'))
        shader_program._set_uniform("camera_pos",
                self.camera.position)

    def _compute_geometry(self, block, bbox_vertices):
        move = get_translate_matrix(*block.LeftEdge)
        dds = (block.RightEdge - block.LeftEdge)
        scale = get_scale_matrix(*dds)

        transformed_box = bbox_vertices.dot(scale.T).dot(move.T).astype("float32")
        return transformed_box

    def _load_textures(self):
        print "Loading textures."
        if len(self.gl_texture_names) == 0:
            self.gl_texture_names = GL.glGenTextures(len(self.blocks))
            if len(self.blocks) == 1:
                self.gl_texture_names = [self.gl_texture_names]
            GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER,
                    GL.GL_LINEAR)

        else:
            for texture in self.gl_texture_names:
                GL.glDeleteTextures([texture])

            self.gl_texture_names = GL.glGenTextures(len(self.blocks))
            if len(self.blocks) == 1:
                self.gl_texture_names = [self.gl_texture_names]

        for block_id in sorted(self.blocks):
            tex_i, block = self.blocks[block_id]
            texture_name = self.gl_texture_names[tex_i]
            dx, dy, dz = block.my_data[0].shape
            n_data = block.my_data[0].copy(order="F").astype("float32")
            n_data = (n_data - self.min_val) / ((self.max_val - self.min_val) * self.diagonal)
            GL.glBindTexture(GL.GL_TEXTURE_3D, texture_name)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_EDGE)
            GL.glTexStorage3D(GL.GL_TEXTURE_3D, 1, GL.GL_R32F,
                *block.my_data[0].shape)
            GL.glTexSubImage3D(GL.GL_TEXTURE_3D, 0, 0, 0, 0, dx, dy, dz,
                        GL.GL_RED, GL.GL_FLOAT, n_data.T)
            GL.glGenerateMipmap(GL.GL_TEXTURE_3D)

class SceneGraph(SceneComponent):
    def __init__(self):
        super(SceneGraph, self).__init__()
        self.collections = []
        self.fbo = None
        self.fb_texture = None
        self.cmap_texture = None
        self.camera = None
        self.shader_program = None
        self.fb_shader_program = None
        self.min_val, self.max_val = 1e60, -1e60
        self.diagonal = 0.0
        self.data_logged = True

        ox, oy, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        self.width = width
        self.height = height

        self.set_shader("passthrough.v")
        self.set_shader("apply_colormap.f")

        self._init_framebuffer()

    def _init_framebuffer(self):
        self._initialize_vertex_array("fb_vbo")
        GL.glBindVertexArray(self.vert_arrays["fb_vbo"])

        quad_attrib = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, quad_attrib)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, FULLSCREEN_QUAD.nbytes,
                        FULLSCREEN_QUAD, GL.GL_STATIC_DRAW)
        GL.glVertexAttribPointer(0, 3, GL.GL_FLOAT, GL.GL_FALSE, 0, None)

        # unbind
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        GL.glBindVertexArray(0)

        self.setup_fb(self.width, self.height)

    def setup_cmap_tex(self):
        '''Creates 1D texture that will hold colormap in framebuffer'''
        self.cmap_texture = GL.glGenTextures(1)   # create target texture
        GL.glBindTexture(GL.GL_TEXTURE_1D, self.cmap_texture)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexParameterf(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameteri(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexImage1D(GL.GL_TEXTURE_1D, 0, GL.GL_RGBA, 256,
                        0, GL.GL_RGBA, GL.GL_FLOAT, self.camera.cmap)
        GL.glBindTexture(GL.GL_TEXTURE_1D, 0)


    def update_cmap_tex(self):
        '''Updates 1D texture with colormap that's used in framebuffer'''
        if self.camera is None or not self.camera.cmap_new:
            return

        if self.cmap_texture is None:
            self.setup_cmap_tex()

        GL.glBindTexture(GL.GL_TEXTURE_1D, self.cmap_texture)
        GL.glTexSubImage1D(GL.GL_TEXTURE_1D, 0, 0, 256,
                           GL.GL_RGBA, GL.GL_FLOAT, self.camera.cmap)
        GL.glBindTexture(GL.GL_TEXTURE_1D, 0)
        self.camera.cmap_new = False


    def setup_fb(self, width, height):
        '''Setups FrameBuffer that will be used as container
           for 1 pass of rendering'''
        # Clean up old FB and Texture
        if self.fb_texture is not None and \
            GL.glIsTexture(self.fb_texture):
                GL.glDeleteTextures([self.fb_texture])
        if self.fbo is not None and GL.glIsFramebuffer(self.fbo):
            GL.glDeleteFramebuffers(1, [self.fbo])


        # initialize FrameBuffer
        self.fbo = GL.glGenFramebuffers(1)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo)

        depthbuffer = GL.glGenRenderbuffers(1)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, depthbuffer)
        GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT32F,
                                 width, height)
        GL.glFramebufferRenderbuffer(
            GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT, GL.GL_RENDERBUFFER,
            depthbuffer
        )
        # end of FrameBuffer initialization

        # generate the texture we render to, and set parameters
        self.fb_texture = GL.glGenTextures(1)   # create target texture
        # bind to new texture, all future texture functions will modify this
        # particular one
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.fb_texture)
        # set how our texture behaves on x,y boundaries
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
        # set how our texture is filtered
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)

        # occupy width x height texture memory, (None at the end == empty
        # image)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA32F, width,
                        height, 0, GL.GL_RGBA, GL.GL_FLOAT, None)

        # --- end texture init

        # Set "fb_texture" as our colour attachement #0
        GL.glFramebufferTexture2D(
            GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D,
            self.fb_texture,
            0 # mipmap level, normally 0
        )

        # verify that everything went well
        status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
        assert status == GL.GL_FRAMEBUFFER_COMPLETE, status


    def add_collection(self, collection):
        r"""Adds a block collection to the scene. Collections must not overlap.

        Although it hasn't been tested, this should allow one to add multiple
        datasets to a single scene.

        Parameters
        ----------
        collection : BlockCollection
            A block collection object representing a data-set. Behavior is
            undefined if any collection overlaps with another.

        """
        self.collections.append(collection)
        self.update_minmax()

    def update_minmax(self):
        self.min_val, self.max_val, self.diagonal = 1e60, -1e60, -1e60
        self.data_logged = False

        for collection in self.collections:
            self.min_val = min(self.min_val, collection.min_val)
            self.max_val = max(self.max_val, collection.max_val)
            # doesn't make sense for multiple collections
            self.diagonal = max(self.diagonal, collection.diagonal)
            self.data_logged = self.data_logged or collection.data_logged

        if self.camera is not None:
            self.camera.update_cmap_minmax(self.min_val, self.max_val,
                                           self.data_logged)

    def set_camera(self, camera):
        r""" Sets the camera orientation for the entire scene.

        Upon calling this function a kd-tree is constructed for each collection.
        This function simply calls the BlockCollection set_camera function on
        each collection in the scene.

        Parameters
        ----------
        camera : Camera

        """
        self.camera = camera
        for collection in self.collections:
            collection.set_camera(camera)

    def _retrieve_framebuffer(self):
        ox, oy, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        debug_buffer = GL.glReadPixels(0, 0, width, height, GL.GL_RGB,
                                       GL.GL_UNSIGNED_BYTE)
        arr = np.fromstring(debug_buffer, "uint8", count = width*height*3)
        return arr.reshape((width, height, 3))

    def run_program(self):
        """ Renders one frame of the scene.

        Renders the scene using the current collection and camera set by calls
        to add_collection and set_camera respectively. Also uses the last shader
        provided to the add_shader_from_file function.

        """

        # get size of current viewport
        ox, oy, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        if (width, height) != (self.width, self.height):
            # size of viewport changed => fb needs to be recreated
            self.setup_fb(width, height)
            self.width = width
            self.width = height

        # Handle colormap
        self.update_cmap_tex()

        # bind to fb
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo)
        # clear the color and depth buffer
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        # render collections to fb
        for collection in self.collections:
            collection.run_program()
        # unbind FB
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)

        # 2 pass of rendering
        GL.glUseProgram(self.program.program)
        GL.glActiveTexture(GL.GL_TEXTURE0)
        # bind to the result of 1 pass
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.fb_texture)
        # Set our "fb_texture" sampler to user Texture Unit 0
        self.program._set_uniform("fb_texture", 0)

        GL.glActiveTexture(GL.GL_TEXTURE1)
        GL.glBindTexture(GL.GL_TEXTURE_1D, self.cmap_texture)
        self.program._set_uniform("cmap", 1)

        scale = (self.max_val - self.min_val) * self.diagonal
        self.program._set_uniform("min_val", self.min_val)
        self.program._set_uniform("scale", scale)
        self.program._set_uniform("cmap_min", self.camera.cmap_min)
        self.program._set_uniform("cmap_max", self.camera.cmap_max)
        if self.data_logged:
            self.program._set_uniform("cmap_log", float(False))
        else:
            self.program._set_uniform("cmap_log", float(self.camera.cmap_log))
        # clear the color and depth buffer
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        # Bind to Vertex array that contains simple quad filling fullscreen,
        # that was defined in __init__()
        GL.glBindVertexArray(self.vert_arrays["fb_vbo"])
        GL.glEnableVertexAttribArray(0)
        # Draw our 2 triangles
        GL.glDrawArrays(GL.GL_TRIANGLES, 0, 6)
        # Clean up
        GL.glDisableVertexAttribArray(0)
        GL.glBindVertexArray(0)
        GL.glBindTexture(GL.GL_TEXTURE_1D, 0)
        GL.glBindTexture(GL.GL_TEXTURE_2D, 0)

    render = run_program
