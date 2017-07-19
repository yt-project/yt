# encoding: utf-8
"""
Interactive Data Visualization classes for Scene, Camera and BlockCollection

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# This is a part of the experimental Interactive Data Visualization

import OpenGL.GL as GL
import OpenGL.GLUT.freeglut as GLUT
from collections import OrderedDict, namedtuple
import matplotlib.cm as cm
import matplotlib.font_manager
from matplotlib.ft2font import FT2Font, LOAD_FORCE_AUTOHINT, LOAD_NO_HINTING, \
     LOAD_DEFAULT, LOAD_NO_AUTOHINT
import numpy as np
import ctypes
import time
import traitlets
import string

from yt import write_bitmap
from yt.config import \
    ytcfg
from yt.utilities.math_utils import \
    get_translate_matrix, \
    get_scale_matrix, \
    get_lookat_matrix, \
    get_perspective_matrix, \
    quaternion_mult, \
    quaternion_to_rotation_matrix, \
    rotation_matrix_to_quaternion
from yt.data_objects.data_containers import \
    YTDataContainer
from yt.data_objects.static_output import \
     Dataset
from yt.utilities.lib.mesh_triangulation import triangulate_mesh
from yt.extern.six import unichr
from .shader_objects import \
    known_shaders, ShaderProgram, ShaderTrait
from .opengl_support import \
    Texture, Texture1D, Texture2D, Texture3D, \
    VertexArray, VertexAttribute, ColormapTexture, \
    Framebuffer

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

class IDVCamera(object):
    '''Camera object used in the Interactive Data Visualization

    Parameters
    ----------

    position : :obj:`!iterable`, or 3 element array in code_length
        The initial position of the camera.
    focus : :obj:`!iterable`, or 3 element array in code_length
        A point in space that the camera is looking at.
    up : :obj:`!iterable`, or 3 element array in code_length
        The 'up' direction for the camera.
    fov : float, optional
        An angle defining field of view in degrees.
    near_plane : float, optional
        The distance to the near plane of the perspective camera.
    far_plane : float, optional
        The distance to the far plane of the perspective camera.
    aspect_ratio: float, optional
        The ratio between the height and the width of the camera's fov.

    '''
    def __init__(self,
                 position=(0.0, 0.0, 1.0),
                 focus=(0.0, 0.0, 0.0),
                 up=(0.0, 1.0, 0.0),
                 fov=45.0, near_plane=0.01, far_plane=20.0,
                 aspect_ratio=8.0/6.0):
        self.position = np.array(position)
        self.focus = np.array(focus)
        self.up = np.array(up)
        self.fov = fov
        self.near_plane = near_plane
        self.far_plane = far_plane
        self.aspect_ratio = aspect_ratio
        
        # set cmap
        cmap = cm.get_cmap(ytcfg.get("yt", "default_colormap"))
        self.cmap = np.array(cmap(np.linspace(0, 1, 256)), dtype=np.float32)
        self.cmap_min = 1e55
        self.cmap_max = -1e55
        self.cmap_log = True
        self.cmap_new = True

        self.view_matrix = np.zeros((4, 4), dtype=np.float32)
        self.projection_matrix = np.zeros((4, 4), dtype=np.float32)
        self.orientation = np.zeros((4, 4), dtype=np.float32)
        self.proj_func = get_perspective_matrix

    def compute_matrices(self):
        '''Regenerate all position, view and projection matrices of the camera.'''
        pass

    def update_orientation(self, start_x, start_y, end_x, end_y):
        '''Change camera orientation matrix using delta of mouse's cursor position
        
        Parameters
        ----------

        start_x : float
            initial cursor position in x direction
        start_y : float
            initial cursor position in y direction
        end_x : float
            final cursor position in x direction
        end_y : float
            final cursor position in y direction

        '''
        pass

    def get_viewpoint(self):
        return self.position

    def get_view_matrix(self):
        return self.view_matrix

    def get_projection_matrix(self):
        return self.projection_matrix

    def update_cmap_minmax(self, minval, maxval, iflog):
        '''Update camera's colormap bounds

        Parameters
        ----------
        
        minval: float
            min color limit used for image scaling
        maxval: float
            max color limit used for image scaling
        iflog: boolean
            Set to True if colormap is using log scale, False for linear scale.
        '''
        self.cmap_log = iflog
        self.cmap_min = minval
        self.cmap_max = maxval

class TrackballCamera(IDVCamera):
    """

    This class implements a basic "Trackball" or "Arcball" camera control system
    that allows for unconstrained 3D rotations without suffering from Gimbal lock.
    Following Ken Shoemake's orginal C implementation (Graphics Gems IV, III.1)
    we project mouse movements onto the unit sphere and use quaternions to
    represent the corresponding rotation.

    See also:
    https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Arcball

    """

    def __init__(self, position=(0.0, 0.0, 1.0), focus=(0.0, 0.0, 0.0),
                 up=(0.0, 1.0, 0.0), fov=45.0, near_plane=0.01, far_plane=20.0,
                 aspect_ratio=8.0/6.0):

        super(TrackballCamera, self).__init__(position=position, focus=focus,
                                              up=up, fov=fov, 
                                              near_plane=near_plane,
                                              far_plane=far_plane,
                                              aspect_ratio=aspect_ratio)

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

class SceneData(traitlets.HasTraits):
    """A class that defines a collection of GPU-managed data.

    This class contains the largest common set of features that can be used 
    OpenGL rendering: a set of vertices and a set of vertex attributes.  Note
    that this is distinct from the shader, which can be swapped out and
    provided with these items.

    """
    name = None
    vertex_array = traitlets.Instance(VertexArray)
    textures = traitlets.List(trait = traitlets.Instance(Texture))

class SceneComponent(traitlets.HasTraits):
    data = traitlets.Instance(SceneData)
    base_quad = traitlets.Instance(SceneData)
    fragment_shader = ShaderTrait(allow_none = True)
    vertex_shader = ShaderTrait(allow_none = True)
    fb = traitlets.Instance(Framebuffer)
    colormap_fragment = ShaderTrait(allow_none = True)
    colormap_vertex = ShaderTrait(allow_none = True)
    colormap = traitlets.Instance(ColormapTexture)
    _program1 = traitlets.Instance(ShaderProgram, allow_none = True)
    _program2 = traitlets.Instance(ShaderProgram, allow_none = True)
    _program1_invalid = True
    _program2_invalid = True

    # These attributes are 
    min_val = traitlets.CFloat(0.0)
    cmin = traitlets.CFloat(0.0)
    cmax = traitlets.CFloat(1.0)
    cmap_log = traitlets.Bool(False)
    scale = traitlets.CFloat(1.0)

    @traitlets.default("fb")
    def _fb_default(self):
        return Framebuffer()

    @traitlets.observe("fragment_shader")
    def _change_fragment(self, change):
        # Even if old/new are the same
        self._program1_invalid = True

    @traitlets.observe("vertex_shader")
    def _change_vertex(self, change):
        # Even if old/new are the same
        self._program1_invalid = True

    @traitlets.observe("colormap_vertex")
    def _change_colormap_vertex(self, change):
        # Even if old/new are the same
        self._program2_invalid = True

    @traitlets.observe("colormap_shader")
    def _change_colormap_fragment(self, change):
        # Even if old/new are the same
        self._program2_invalid = True

    @traitlets.default("colormap")
    def _default_colormap(self):
        cm = ColormapTexture()
        cm.colormap_name = "arbre"
        return cm

    @traitlets.default("base_quad")
    def _default_base_quad(self):
        bq = SceneData(name = "fullscreen_quad", 
                       vertex_array = VertexArray(name = "tri", each = 6),
        )
        fq = FULLSCREEN_QUAD.reshape((6, 3), order="C")
        bq.vertex_array.attributes.append(VertexAttribute(
            name = "vertexPosition_modelspace", data = fq
        ))
        return bq

    @property
    def program1(self):
        if self._program1_invalid:
            if self._program1 is not None:
                self._program1.delete_program()
            self._program1 = ShaderProgram(self.vertex_shader,
                self.fragment_shader)
            self._program1_invalid = False
        return self._program1

    @property
    def program2(self):
        if self._program2_invalid:
            if self._program2 is not None:
                self._program2.delete_program()
            # The vertex shader will always be the same.
            # The fragment shader will change based on whether we are
            # colormapping or not.
            self._program2 = ShaderProgram(self.colormap_vertex,
                    self.colormap_fragment)
            self._program2_invalid = False
        return self._program2

    def run_program(self, scene):
        self.init_draw(scene)
        with self.fb.bind(True):
            with self.program1.enable() as p:
                self._set_uniforms(scene, p)
                with self.data.vertex_array.bind(p):
                    self.draw(scene, p)
        with self.colormap.bind(0):
            with self.fb.input_bind(1, 2):
                with self.program2.enable() as p2:
                    self.init_fb_draw(scene)
                    p2._set_uniform("cmap", 0)
                    p2._set_uniform("fb_texture", 1)
                    p2._set_uniform("db_texture", 2)
                    p2._set_uniform("min_val", self.min_val)
                    p2._set_uniform("scale", self.scale)
                    p2._set_uniform("cmap_min", self.cmin)
                    p2._set_uniform("cmap_max", self.cmax)
                    p2._set_uniform("cmap_log", float(self.cmap_log))
                    with self.base_quad.vertex_array.bind(p2):
                        GL.glDrawArrays(GL.GL_TRIANGLES, 0, 6)
                
    def draw(self, scene, program):
        raise NotImplementedError

    def init_draw(self, scene):
        return

    def init_fb_draw(self, scene):
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthMask(GL.GL_TRUE)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        return

class SceneAnnotation(SceneComponent):
    pass

# This is drawn in part from
#  https://learnopengl.com/#!In-Practice/Text-Rendering
Character = namedtuple('Character',
        ['texture', 'vbo_offset', 'hori_advance', 'vert_advance']
)

class FontTrait(traitlets.TraitType):
    info_text = "A font instance from matplotlib"

    def validate(self, obj, value):
        if isinstance(value, str):
            try:
                font_fn = matplotlib.font_manager.findfont(value)
                value = matplotlib.font_manager.get_font(font_fn)
            except FileNotFoundError:
                self.error(obj, value)
        return value

class TextCharacters(SceneData):
    characters = traitlets.Dict(trait=traitlets.Instance(Character))
    name = "text_overlay"
    font = FontTrait("DejaVu Sans")
    font_size = traitlets.CInt(32)

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        return VertexArray(name = "char_info", each = 6)

    def build_textures(self):
        # This doesn't check if the textures have already been built
        self.font.set_size(self.font_size, 200)
        chars = [ord(_) for _ in string.printable]
        tex_ids = GL.glGenTextures(len(chars))
        vert = []
        for i, (tex_id, char_code) in enumerate(zip(tex_ids, chars)):
            self.font.clear()
            self.font.set_text(unichr(char_code),
                    flags = LOAD_FORCE_AUTOHINT)
            self.font.draw_glyphs_to_bitmap(antialiased = True)
            glyph = self.font.load_char(char_code)
            x0, y0, x1, y1 = glyph.bbox
            bitmap = self.font.get_image().astype(">f4")/255.0
            dx = 1.0/bitmap.shape[0]
            dy = 1.0/bitmap.shape[1]
            triangles = np.array([[x0, y1, 0.0 + dx/2.0, 0.0 + dy/2.0],
                                  [x0, y0, 0.0 + dx/2.0, 1.0 - dy/2.0],
                                  [x1, y0, 1.0 - dx/2.0, 1.0 - dy/2.0],
                                  [x0, y1, 0.0 + dx/2.0, 0.0 + dy/2.0],
                                  [x1, y0, 1.0 - dx/2.0, 1.0 - dy/2.0],
                                  [x1, y1, 1.0 - dx/2.0, 0.0 + dy/2.0]],
                                  dtype="<f4")
            vert.append(triangles)
            texture = Texture2D(texture_name = tex_id,
                                data = bitmap,
                                boundary_x = "clamp", 
                                boundary_y = "clamp")
            # I can't find information as to why horiAdvance is a
            # factor of 8 larger than the other factors.  I assume it
            # is referenced somewhere, but I cannot find it.
            self.characters[unichr(char_code)] = Character(texture,
                    i, glyph.horiAdvance/8., glyph.vertAdvance)
        vert = np.concatenate(vert)
        self.vertex_array.attributes.append(VertexAttribute(
            name = "quad_vertex", data = vert.astype("<f4")))

class TextAnnotation(SceneAnnotation):

    data = traitlets.Instance(TextCharacters)
    text = traitlets.CUnicode()
    draw_instructions = traitlets.List()
    origin = traitlets.Tuple(traitlets.CFloat(), traitlets.CFloat(),
            default_value = (-1, -1))
    scale = traitlets.CFloat(1.0)

    @traitlets.observe("text")
    def _observe_text(self, change):
        text = change['new']
        lines = text.split("\n")
        draw_instructions = []
        y = 0
        for line in reversed(lines):
            x = 0
            dy = 0
            for c in line:
                e = self.data.characters[c]
                draw_instructions.append((x, y,
                    e.texture, e.vbo_offset))
                dy = max(dy, e.vert_advance)
                x += e.hori_advance
            y += dy
        self.draw_instructions = draw_instructions

    def _set_uniforms(self, scene, shader_program):
        pass

    def draw(self, scene, program):
        viewport = np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype="f4")
        program._set_uniform("viewport", viewport)
        each = self.data.vertex_array.each
        for x, y, tex, vbo_offset in self.draw_instructions:
            with tex.bind(0):
                program._set_uniform("x_offset", float(x))
                program._set_uniform("y_offset", float(y))
                program._set_uniform("x_origin", self.origin[0])
                program._set_uniform("y_origin", self.origin[1])
                program._set_uniform("scale", self.scale)
                GL.glDrawArrays(GL.GL_TRIANGLES, vbo_offset*each, each)

    def init_draw(self, scene):
        GL.glDisable(GL.GL_BLEND)

    def init_fb_draw(self, scene):
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

class BlockCollection(SceneData):
    name = "block_collection"
    data_source = traitlets.Instance(YTDataContainer)
    texture_objects = traitlets.List(trait = traitlets.Instance(Texture3D))
    blocks = traitlets.Dict(default_value = ())
    scale = traitlets.Bool(False)

    def _init_blending(self):
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendColor(1.0, 1.0, 1.0, 1.0)
        GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
        GL.glBlendEquation(GL.GL_MAX)

    def set_fields_log(self, log_field):
        """Switch between a logarithmic and a linear scale for the data.

        Parameters
        ----------

        log_field : boolean
            If set to True log10 will be applied to data before passing it to GPU.

        """
        self.add_data(self.data_source, self.data_source.tiles.fields[0], log_field)

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        return VertexArray(name = "block_info", each = 36)

    def add_data(self, field, log_field=True):
        r"""Adds a source of data for the block collection.

        Given a `data_source` and a `field` to populate from, adds the data
        to the block collection so that is able to be rendered.

        Parameters
        ----------
        data_source : YTRegion
            A YTRegion object to use as a data source.
        field : string
            A field to populate from.
        log_field : boolean, optional
            If set to True log10 will be applied to data before passing it to GPU.
        """
        self.data_source.tiles.set_fields([field], [log_field], no_ghost=False)
        self.data_logged = log_field
        # Every time we change our data source, we wipe all existing ones.
        # We now set up our vertices into our current data source.
        vert, dx, le, re = [], [], [], []
        self.min_val = 1e60
        self.max_val = -1e60
        if self.scale:
            left_min = np.ones(3, "f8") * np.inf
            right_max = np.ones(3, "f8") * -np.inf
            for block in self.data_source.tiles.traverse():
                np.minimum(left_min, block.LeftEdge, left_min)
                np.maximum(right_max, block.LeftEdge, right_max)
            scale = right_max.max() - left_min.min()
            for block in self.data_source.tiles.traverse():
                block.LeftEdge -= left_min
                block.LeftEdge /= scale
                block.RightEdge -= left_min
                block.RightEdge /= scale
        for i, block in enumerate(self.data_source.tiles.traverse()):
            self.min_val = min(self.min_val, np.nanmin(block.my_data[0].min()))
            self.max_val = max(self.max_val, np.nanmax(block.my_data[0].max()))
            self.blocks[id(block)] = (i, block)
            vert.append(self._compute_geometry(block, bbox_vertices))
            dds = (block.RightEdge - block.LeftEdge)/block.my_data[0].shape
            n = int(vert[-1].size) // 4
            dx.append([dds.astype('f4') for _ in range(n)])
            le.append([block.LeftEdge.astype('f4') for _ in range(n)])
            re.append([block.RightEdge.astype('f4') for _ in range(n)])

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

        self.vertex_array.attributes.append(VertexAttribute(
            name = "model_vertex", data = vert))
        self.vertex_array.attributes.append(VertexAttribute(
            name = "in_dx", data = dx))
        self.vertex_array.attributes.append(VertexAttribute(
            name = "in_left_edge", data = le))
        self.vertex_array.attributes.append(VertexAttribute(
            name = "in_right_edge", data = re))

        # Now we set up our textures
        self._load_textures()

    def viewpoint_iter(self, camera):
        for block in self.data_source.tiles.traverse(viewpoint = camera.position):
            tex_i, _ = self.blocks[id(block)]
            yield tex_i, self.texture_objects[tex_i]

    def _compute_geometry(self, block, bbox_vertices):
        move = get_translate_matrix(*block.LeftEdge)
        dds = (block.RightEdge - block.LeftEdge)
        scale = get_scale_matrix(*dds)

        transformed_box = bbox_vertices.dot(scale.T).dot(move.T).astype("float32")
        return transformed_box

    def _load_textures(self):
        for block_id in sorted(self.blocks):
            tex_i, block = self.blocks[block_id]
            n_data = block.my_data[0].copy(order="F").astype("float32")
            n_data = (n_data - self.min_val) / ((self.max_val - self.min_val) * self.diagonal)
            tex = Texture3D(data = n_data)
            self.texture_objects.append(tex)

class BlockRendering(SceneComponent):
    '''
    A class that renders block data.  It may do this in one of several ways,
    including mesh outline.  This allows us to render a single collection of
    blocks multiple times in a single scene and to separate out the memory
    handling from the display.
    '''
    data = traitlets.Instance(BlockCollection)

    def draw(self, scene, program):
        each = self.data.vertex_array.each
        for tex_ind, texture in self.data.viewpoint_iter(scene.camera):
            with texture.bind(target = 0):
                GL.glDrawArrays(GL.GL_TRIANGLES, tex_ind*each, each)

    def _set_uniforms(self, scene, shader_program):
        cam = scene.camera
        shader_program._set_uniform("projection",
                cam.get_projection_matrix())
        shader_program._set_uniform("modelview",
                cam.get_view_matrix())
        shader_program._set_uniform("viewport",
                np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype = 'f4'))
        shader_program._set_uniform("camera_pos",
                cam.position)
        shader_program._set_uniform("box_width", 1.0)

    def init_draw(self, scene):
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthFunc(GL.GL_LESS)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)

class MeshSceneComponent(object):
    '''

    A scene component for representing unstructured mesh data.

    '''

    def __init__(self, data_source, field):
        super(MeshSceneComponent, self).__init__()
        self.set_shader("mesh.v")
        self.set_shader("mesh.f")

        self.data_source = None
        self.redraw = True

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthFunc(GL.GL_LESS)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)

        vertices, data, indices = self.get_mesh_data(data_source, field)

        self._initialize_vertex_array("mesh_info")
        GL.glBindVertexArray(self.vert_arrays["mesh_info"])

        self.add_vert_attrib("vertex_buffer", vertices, vertices.size)
        self.add_vert_attrib("data_buffer", data, data.size)

        self.vert_attrib["element_buffer"] = (GL.glGenBuffers(1), indices.size)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.vert_attrib["element_buffer"][0])
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_STATIC_DRAW)

        self.transform_matrix = GL.glGetUniformLocation(self.program.program,
                                                        "model_to_clip")

        self.cmin = data.min()
        self.cmax = data.max()

    def set_camera(self, camera):
        r""" Sets the camera orientation for the entire scene.

        Parameters
        ----------
        camera : :class:`~yt.visualization.volume_rendering.interactive_vr.IDVCamera`

        """
        self.camera = camera
        self.camera.cmap_min = float(self.cmin)
        self.camera.cmap_max = float(self.cmax)
        self.redraw = True

    def get_mesh_data(self, data_source, field):
        """
        
        This reads the mesh data into a form that can be fed in to OpenGL.
        
        """

        # get mesh information
        try:
            ftype, fname = field
            mesh_id = int(ftype[-1])
        except ValueError:
            mesh_id = 0

        mesh = data_source.ds.index.meshes[mesh_id-1]
        offset = mesh._index_offset
        vertices = mesh.connectivity_coords
        indices  = mesh.connectivity_indices - offset

        data = data_source[field]

        return triangulate_mesh(vertices, data, indices)

    def run_program(self):
        """ Renders one frame of the scene. """
        with self.program.enable():

            # Handle colormap
            self.update_cmap_tex()

            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
            projection_matrix = self.camera.projection_matrix
            view_matrix = self.camera.view_matrix
            model_to_clip = np.dot(projection_matrix, view_matrix)
            GL.glUniformMatrix4fv(self.transform_matrix, 1, True, model_to_clip)

            GL.glActiveTexture(GL.GL_TEXTURE1)
            GL.glBindTexture(GL.GL_TEXTURE_1D, self.cmap_texture)

            self.program._set_uniform("cmap", 0)
            self.program._set_uniform("cmap_min", self.camera.cmap_min)
            self.program._set_uniform("cmap_max", self.camera.cmap_max)

            GL.glEnableVertexAttribArray(0)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vert_attrib["vertex_buffer"][0])
            GL.glVertexAttribPointer(0, 3, GL.GL_FLOAT, False, 0, ctypes.c_void_p(0))

            GL.glEnableVertexAttribArray(1)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vert_attrib["data_buffer"][0])
            GL.glVertexAttribPointer(1, 1, GL.GL_FLOAT, False, 0, ctypes.c_void_p(0))

            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.vert_attrib["element_buffer"][0])
            GL.glDrawElements(GL.GL_TRIANGLES, self.vert_attrib["element_buffer"][1],
                              GL.GL_UNSIGNED_INT, ctypes.c_void_p(0))

            GL.glDisableVertexAttribArray(0)
            GL.glDisableVertexAttribArray(1)

    render = run_program

class SceneGraph(traitlets.HasTraits):
    components = traitlets.List(trait = traitlets.Instance(SceneComponent),
            default_value = [])
    annotations = traitlets.List(trait = traitlets.Instance(SceneAnnotation),
            default_value = [])
    data_objects = traitlets.List(trait = traitlets.Instance(SceneData),
            default_value = [])
    camera = traitlets.Instance(IDVCamera)
    ds = traitlets.Instance(Dataset)

    def render(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for component in self.components:
            component.run_program(self)
        for annotation in self.annotations:
            annotation.run_program(self)

    def set_camera(self, camera):
        self.camera = camera

    def update_minmax(self):
        pass
