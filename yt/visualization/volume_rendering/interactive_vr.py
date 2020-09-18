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

import contextlib
import string
from collections import defaultdict, namedtuple
from math import ceil, floor

import matplotlib.font_manager
import numpy as np
import traitlets
import traittypes
from matplotlib.ft2font import LOAD_FORCE_AUTOHINT
from OpenGL import GL

from yt.data_objects.api import Dataset
from yt.data_objects.data_containers import YTDataContainer
from yt.utilities.lib.mesh_triangulation import triangulate_mesh
from yt.utilities.lib.misc_utilities import update_orientation
from yt.utilities.math_utils import (
    get_lookat_matrix,
    get_perspective_matrix,
    get_scale_matrix,
    get_translate_matrix,
    quaternion_mult,
    quaternion_to_rotation_matrix,
    rotation_matrix_to_quaternion,
)
from yt.utilities.traitlets_support import YTPositionTrait, ndarray_ro, ndarray_shape

from .input_events import EventCollection
from .opengl_support import (
    ColormapTexture,
    Framebuffer,
    Texture,
    Texture1D,
    Texture2D,
    Texture3D,
    TransferFunctionTexture,
    VertexArray,
    VertexAttribute,
)
from .shader_objects import (
    ShaderProgram,
    ShaderTrait,
    component_shaders,
    default_shader_combos,
)

bbox_vertices = np.array(
    [
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0],
    ],
    dtype=np.float32,
)

FULLSCREEN_QUAD = np.array(
    [
        -1.0,
        -1.0,
        0.0,
        +1.0,
        -1.0,
        0.0,
        -1.0,
        +1.0,
        0.0,
        -1.0,
        +1.0,
        0.0,
        +1.0,
        -1.0,
        0.0,
        +1.0,
        +1.0,
        0.0,
    ],
    dtype=np.float32,
)


def compute_box_geometry(left_edge, right_edge):
    move = get_translate_matrix(*left_edge)
    width = right_edge - left_edge
    scale = get_scale_matrix(*width)

    transformed_box = bbox_vertices.dot(scale.T).dot(move.T).astype("float32")
    return transformed_box


class IDVCamera(traitlets.HasTraits):
    """Camera object used in the Interactive Data Visualization

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

    """

    # We have to be careful about some of these, as it's possible in-place
    # operations won't trigger our observation.
    position = YTPositionTrait([0.0, 0.0, 1.0])
    focus = YTPositionTrait([0.0, 0.0, 0.0])
    up = traittypes.Array(np.array([0.0, 0.0, 1.0])).valid(
        ndarray_shape(3), ndarray_ro()
    )
    fov = traitlets.Float(45.0)
    near_plane = traitlets.Float(0.01)
    far_plane = traitlets.Float(20.0)
    aspect_ratio = traitlets.Float(8.0 / 6.0)

    projection_matrix = traittypes.Array(np.zeros((4, 4))).valid(
        ndarray_shape(4, 4), ndarray_ro()
    )
    view_matrix = traittypes.Array(np.zeros((4, 4))).valid(
        ndarray_shape(4, 4), ndarray_ro()
    )
    orientation = traittypes.Array(np.zeros(4)).valid(ndarray_shape(4), ndarray_ro())

    held = traitlets.Bool(False)

    @contextlib.contextmanager
    def hold_traits(self, func):
        # for some reason, hold_trait_notifications doesn't seem to work here.
        # So, we use this to block.  We also do not want to pass the
        # notifications once completed.
        if not self.held:
            self.held = True
            func()
            self.held = False
        yield

    @traitlets.default("up")
    def _default_up(self):
        return np.array([0.0, 1.0, 0.0])

    @traitlets.observe(
        "position",
        "focus",
        "up",
        "fov",
        "near_plane",
        "far_plane",
        "aspect_ratio",
        "orientation",
    )
    def compute_matrices(self, change=None):
        """Regenerate all position, view and projection matrices of the camera."""
        with self.hold_traits(self._compute_matrices):
            pass

    def update_orientation(self, start_x, start_y, end_x, end_y):
        """Change camera orientation matrix using delta of mouse's cursor position

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

        """
        pass


class TrackballCamera(IDVCamera):
    """

    This class implements a basic "Trackball" or "Arcball" camera control system
    that allows for unconstrained 3D rotations without suffering from Gimbal lock.
    Following Ken Shoemake's original C implementation (Graphics Gems IV, III.1)
    we project mouse movements onto the unit sphere and use quaternions to
    represent the corresponding rotation.

    See also:
    https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Arcball

    """

    @property
    def proj_func(self):
        return get_perspective_matrix

    @traitlets.default("orientation")
    def _orientation_default(self):
        rotation_matrix = self.view_matrix[0:3, 0:3]
        return rotation_matrix_to_quaternion(rotation_matrix)

    @traitlets.default("view_matrix")
    def _default_view_matrix(self):
        return get_lookat_matrix(self.position, self.focus, self.up)

    def _map_to_surface(self, mouse_x, mouse_y):
        # right now this just maps to the surface of the unit sphere
        x, y = mouse_x, mouse_y
        mag = np.sqrt(x * x + y * y)
        if mag > 1.0:
            x /= mag
            y /= mag
            z = 0.0
        else:
            z = np.sqrt(1.0 - mag ** 2)
        return np.array([x, -y, z])

    def update_orientation(self, start_x, start_y, end_x, end_y):
        if 0:
            old = self._map_to_surface(start_x, start_y)
            new = self._map_to_surface(end_x, end_y)

            # dot product controls the angle of the rotation
            w = old[0] * new[0] + old[1] * new[1] + old[2] * new[2]

            # cross product gives the rotation axis
            x = old[1] * new[2] - old[2] * new[1]
            y = old[2] * new[0] - old[0] * new[2]
            z = old[0] * new[1] - old[1] * new[0]

            q = np.array([w, x, y, z])

            # renormalize to prevent floating point issues
            mag = np.sqrt(w ** 2 + x ** 2 + y ** 2 + z ** 2)
            q /= mag

            _ = quaternion_mult(self.orientation, q)
        self.orientation = update_orientation(
            self.orientation, start_x, start_y, end_x, end_y
        )

        rotation_matrix = quaternion_to_rotation_matrix(self.orientation)
        dp = np.linalg.norm(self.position - self.focus) * rotation_matrix[2]
        self.position = dp + self.focus
        self.up = rotation_matrix[1]

        self.view_matrix = get_lookat_matrix(self.position, self.focus, self.up)

        self.projection_matrix = self.proj_func(
            self.fov, self.aspect_ratio, self.near_plane, self.far_plane
        )

    def _update_matrices(self):

        self.view_matrix = get_lookat_matrix(self.position, self.focus, self.up)

        self.projection_matrix = self.proj_func(
            self.fov, self.aspect_ratio, self.near_plane, self.far_plane
        )

    def offsetPosition(self, dPos=None):
        if dPos is None:
            dPos = np.array([0.0, 0.0, 0.0])
        self.position += dPos
        self.view_matrix = get_lookat_matrix(self.position, self.focus, self.up)

    def _compute_matrices(self):
        pass


class SceneData(traitlets.HasTraits):
    """A class that defines a collection of GPU-managed data.

    This class contains the largest common set of features that can be used
    OpenGL rendering: a set of vertices and a set of vertex attributes.  Note
    that this is distinct from the shader, which can be swapped out and
    provided with these items.

    """

    name = None
    vertex_array = traitlets.Instance(VertexArray)
    textures = traitlets.List(trait=traitlets.Instance(Texture))

    min_val = traitlets.CFloat(0.0)
    max_val = traitlets.CFloat(1.0)


class SceneComponent(traitlets.HasTraits):
    data = traitlets.Instance(SceneData)
    base_quad = traitlets.Instance(SceneData)
    events = traitlets.Instance(EventCollection)
    name = "undefined"
    priority = traitlets.CInt(0)
    visible = traitlets.Bool(True)
    display_bounds = traitlets.Tuple((0.0, 1.0, 0.0, 1.0), trait=traitlets.CFloat())
    clear_region = traitlets.Bool(False)

    render_method = traitlets.Unicode(allow_none=True)
    fragment_shader = ShaderTrait(allow_none=True).tag(shader_type="fragment")
    vertex_shader = ShaderTrait(allow_none=True).tag(shader_type="vertex")
    fb = traitlets.Instance(Framebuffer)
    colormap_fragment = ShaderTrait(allow_none=True).tag(shader_type="fragment")
    colormap_vertex = ShaderTrait(allow_none=True).tag(shader_type="vertex")
    colormap = traitlets.Instance(ColormapTexture)
    _program1 = traitlets.Instance(ShaderProgram, allow_none=True)
    _program2 = traitlets.Instance(ShaderProgram, allow_none=True)
    _program1_invalid = True
    _program2_invalid = True

    # These attributes are
    cmap_min = traitlets.CFloat(None, allow_none=True)
    cmap_max = traitlets.CFloat(None, allow_none=True)
    cmap_log = traitlets.Bool(True)
    scale = traitlets.CFloat(1.0)

    @traitlets.observe("display_bounds")
    def _change_display_bounds(self, change):
        # We need to update the framebuffer if the width or height has changed
        # Same thing is true if the total pixel size has changed, but that is
        # not doable from in here.
        if change["old"] == traitlets.Undefined:
            return
        old_width = change["old"][1] - change["old"][0]
        old_height = change["old"][3] - change["old"][2]
        new_width = change["new"][1] - change["new"][0]
        new_height = change["new"][3] - change["new"][2]
        if old_width != new_width or old_height != new_height:
            self.fb = Framebuffer()

    def render_gui(self, imgui, renderer):
        changed, self.visible = imgui.checkbox("Visible", self.visible)
        return changed

    @traitlets.default("render_method")
    def _default_render_method(self):
        return default_shader_combos[self.name]

    @traitlets.observe("render_method")
    def _change_render_method(self, change):
        new_combo = component_shaders[self.name][change["new"]]
        with self.hold_trait_notifications():
            self.vertex_shader = new_combo["first_vertex"]
            self.fragment_shader = new_combo["first_fragment"]
            self.colormap_vertex = new_combo["second_vertex"]
            self.colormap_fragment = new_combo["second_fragment"]

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

    @traitlets.observe("colormap_fragment")
    def _change_colormap_fragment(self, change):
        # Even if old/new are the same
        self._program2_invalid = True

    @traitlets.default("colormap")
    def _default_colormap(self):
        cm = ColormapTexture()
        cm.colormap_name = "arbre"
        return cm

    @traitlets.default("vertex_shader")
    def _vertex_shader_default(self):
        return component_shaders[self.name][self.render_method]["first_vertex"]

    @traitlets.default("fragment_shader")
    def _fragment_shader_default(self):
        return component_shaders[self.name][self.render_method]["first_fragment"]

    @traitlets.default("colormap_vertex")
    def _colormap_vertex_default(self):
        return component_shaders[self.name][self.render_method]["second_vertex"]

    @traitlets.default("colormap_fragment")
    def _colormap_fragment_default(self):
        return component_shaders[self.name][self.render_method]["second_fragment"]

    @traitlets.default("base_quad")
    def _default_base_quad(self):
        bq = SceneData(
            name="fullscreen_quad", vertex_array=VertexArray(name="tri", each=6),
        )
        fq = FULLSCREEN_QUAD.reshape((6, 3), order="C")
        bq.vertex_array.attributes.append(
            VertexAttribute(name="vertexPosition_modelspace", data=fq)
        )
        return bq

    @property
    def program1(self):
        if self._program1_invalid:
            if self._program1 is not None:
                self._program1.delete_program()
            self._program1 = ShaderProgram(self.vertex_shader, self.fragment_shader)
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
            self._program2 = ShaderProgram(self.colormap_vertex, self.colormap_fragment)
            self._program2_invalid = False
        return self._program2

    def run_program(self, scene):
        # Store this info, because we need to render into a framebuffer that is the
        # right size.
        x0, y0, w, h = GL.glGetIntegerv(GL.GL_VIEWPORT)
        GL.glViewport(0, 0, w, h)
        if not self.visible:
            return
        with self.fb.bind(True):
            with self.program1.enable() as p:
                self._set_uniforms(scene, p)
                with self.data.vertex_array.bind(p):
                    self.draw(scene, p)
        if self.cmap_min is None or self.cmap_max is None:
            data = self.fb.data
            data = data[data[:, :, 3] > 0][:, 0]
            if self.cmap_min is None and data.size > 0:
                self.cmap_min = cmap_min = data.min()
            if self.cmap_max is None and data.size > 0:
                self.cmap_max = cmap_max = data.max()
            if data.size == 0:
                cmap_min = 0.0
                cmap_max = 1.0
        else:
            cmap_min = float(self.cmap_min)
            cmap_max = float(self.cmap_max)
        with self.colormap.bind(0):
            with self.fb.input_bind(1, 2):
                with self.program2.enable() as p2:
                    with scene.bind_buffer():
                        p2._set_uniform("cmap", 0)
                        p2._set_uniform("fb_texture", 1)
                        p2._set_uniform("db_texture", 2)
                        # Note that we use cmap_min/cmap_max, not
                        # self.cmap_min/self.cmap_max.
                        p2._set_uniform("cmap_min", cmap_min)
                        p2._set_uniform("cmap_max", cmap_max)
                        p2._set_uniform("cmap_log", float(self.cmap_log))
                        with self.base_quad.vertex_array.bind(p2):
                            # Now we do our viewport globally, not just within
                            # the framebuffer
                            GL.glViewport(x0, y0, w, h)
                            GL.glDrawArrays(GL.GL_TRIANGLES, 0, 6)

    def draw(self, scene, program):
        raise NotImplementedError


class SceneAnnotation(SceneComponent):
    pass


# This is drawn in part from
#  https://learnopengl.com/#!In-Practice/Text-Rendering
Character = namedtuple(
    "Character", ["texture", "vbo_offset", "hori_advance", "vert_advance"]
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
        return VertexArray(name="char_info", each=6)

    def build_textures(self):
        # This doesn't check if the textures have already been built
        self.font.set_size(self.font_size, 200)
        chars = [ord(_) for _ in string.printable]
        tex_ids = GL.glGenTextures(len(chars))
        vert = []
        for i, (tex_id, char_code) in enumerate(zip(tex_ids, chars)):
            self.font.clear()
            self.font.set_text(chr(char_code), flags=LOAD_FORCE_AUTOHINT)
            self.font.draw_glyphs_to_bitmap(antialiased=True)
            glyph = self.font.load_char(char_code)
            x0, y0, x1, y1 = glyph.bbox
            bitmap = self.font.get_image().astype(">f4") / 255.0
            dx = 1.0 / bitmap.shape[0]
            dy = 1.0 / bitmap.shape[1]
            triangles = np.array(
                [
                    [x0, y1, 0.0 + dx / 2.0, 0.0 + dy / 2.0],
                    [x0, y0, 0.0 + dx / 2.0, 1.0 - dy / 2.0],
                    [x1, y0, 1.0 - dx / 2.0, 1.0 - dy / 2.0],
                    [x0, y1, 0.0 + dx / 2.0, 0.0 + dy / 2.0],
                    [x1, y0, 1.0 - dx / 2.0, 1.0 - dy / 2.0],
                    [x1, y1, 1.0 - dx / 2.0, 0.0 + dy / 2.0],
                ],
                dtype="<f4",
            )
            vert.append(triangles)
            texture = Texture2D(
                texture_name=tex_id, data=bitmap, boundary_x="clamp", boundary_y="clamp"
            )
            # I can't find information as to why horiAdvance is a
            # factor of 8 larger than the other factors.  I assume it
            # is referenced somewhere, but I cannot find it.
            self.characters[chr(char_code)] = Character(
                texture, i, glyph.horiAdvance / 8.0, glyph.vertAdvance
            )
        vert = np.concatenate(vert)
        self.vertex_array.attributes.append(
            VertexAttribute(name="quad_vertex", data=vert.astype("<f4"))
        )


class TextAnnotation(SceneAnnotation):

    name = "text_annotation"
    data = traitlets.Instance(TextCharacters)
    text = traitlets.CUnicode()
    draw_instructions = traitlets.List()
    origin = traitlets.Tuple(
        traitlets.CFloat(), traitlets.CFloat(), default_value=(-1, -1)
    )
    scale = traitlets.CFloat(1.0)

    @traitlets.observe("text")
    def _observe_text(self, change):
        text = change["new"]
        lines = text.split("\n")
        draw_instructions = []
        y = 0
        for line in reversed(lines):
            x = 0
            dy = 0
            for c in line:
                e = self.data.characters[c]
                draw_instructions.append((x, y, e.texture, e.vbo_offset))
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
                GL.glDrawArrays(GL.GL_TRIANGLES, vbo_offset * each, each)


class LineData(SceneData):
    name = "line_data"
    n_values = traitlets.CInt()

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        return VertexArray(name="vertices", each=6)

    def add_data(self, lines):
        assert lines.shape[1] == 4
        x_coord = np.mgrid[0.0 : 1.0 : lines.shape[0] * 1j].astype("f4")
        x_coord = x_coord.reshape((-1, 1))
        self.n_vertices = lines.shape[0]
        self.vertex_array.attributes.append(
            VertexAttribute(name="rgba_values", data=lines)
        )
        self.vertex_array.attributes.append(
            VertexAttribute(name="x_coord", data=x_coord)
        )


class BlockCollection(SceneData):
    name = "block_collection"
    data_source = traitlets.Instance(YTDataContainer)
    texture_objects = traitlets.Dict(trait=traitlets.Instance(Texture3D))
    bitmap_objects = traitlets.Dict(trait=traitlets.Instance(Texture3D))
    blocks = traitlets.Dict(default_value=())
    scale = traitlets.Bool(False)
    blocks_by_grid = traitlets.Instance(defaultdict, (list,))
    grids_by_block = traitlets.Dict(default_value=())

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        return VertexArray(name="block_info", each=36)

    def add_data(self, field, no_ghost=False):
        r"""Adds a source of data for the block collection.

        Given a `data_source` and a `field` to populate from, adds the data
        to the block collection so that is able to be rendered.

        Parameters
        ----------
        data_source : YTRegion
            A YTRegion object to use as a data source.
        field : string
            A field to populate from.
        no_ghost : bool (False)
            Should we speed things up by skipping ghost zone generation?
        """
        self.data_source.tiles.set_fields([field], [False], no_ghost=no_ghost)
        # Every time we change our data source, we wipe all existing ones.
        # We now set up our vertices into our current data source.
        vert, dx, le, re = [], [], [], []
        self.min_val = +np.inf
        self.max_val = -np.inf
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
            self.min_val = min(self.min_val, np.nanmin(np.abs(block.my_data[0]).min()))
            self.max_val = max(self.max_val, np.nanmax(np.abs(block.my_data[0]).max()))
            self.blocks[id(block)] = (i, block)
            vert.append(compute_box_geometry(block.LeftEdge, block.RightEdge))
            dds = (block.RightEdge - block.LeftEdge) / block.my_data[0].shape
            n = int(vert[-1].size) // 4
            dx.append([dds.astype("f4") for _ in range(n)])
            le.append([block.LeftEdge.astype("f4") for _ in range(n)])
            re.append([block.RightEdge.astype("f4") for _ in range(n)])
        for (g, node, (sl, _dims, _gi)) in self.data_source.tiles.slice_traverse():
            block = node.data
            self.blocks_by_grid[g.id - g._id_offset].append((id(block), i))
            self.grids_by_block[id(node.data)] = (g.id - g._id_offset, sl)

        if hasattr(self.min_val, "in_units"):
            self.min_val = self.min_val.d
        if hasattr(self.max_val, "in_units"):
            self.max_val = self.max_val.d

        LE = np.array([b.LeftEdge for i, b in self.blocks.values()]).min(axis=0)
        RE = np.array([b.RightEdge for i, b in self.blocks.values()]).max(axis=0)
        self.diagonal = np.sqrt(((RE - LE) ** 2).sum())
        # Now we set up our buffer
        vert = np.concatenate(vert)
        dx = np.concatenate(dx)
        le = np.concatenate(le)
        re = np.concatenate(re)

        self.vertex_array.attributes.append(
            VertexAttribute(name="model_vertex", data=vert)
        )
        self.vertex_array.attributes.append(VertexAttribute(name="in_dx", data=dx))
        self.vertex_array.attributes.append(
            VertexAttribute(name="in_left_edge", data=le)
        )
        self.vertex_array.attributes.append(
            VertexAttribute(name="in_right_edge", data=re)
        )

        # Now we set up our textures
        self._load_textures()

    def viewpoint_iter(self, camera):
        for block in self.data_source.tiles.traverse(viewpoint=camera.position):
            vbo_i, _ = self.blocks[id(block)]
            yield (vbo_i, self.texture_objects[vbo_i], self.bitmap_objects[vbo_i])

    def filter_callback(self, callback):
        # This is not efficient.  It calls it once for each node in a grid.
        # We do this the slow way because of the problem of ordering the way we
        # iterate over the grids and nodes.  This can be fixed at some point.
        for g_ind in self.blocks_by_grid:
            blocks = self.blocks_by_grid[g_ind]
            # Does this need an offset?
            grid = self.data_source.index.grids[g_ind]
            new_bitmap = callback(grid).astype("uint8")
            for b_id, _ in blocks:
                _, sl = self.grids_by_block[b_id]
                vbo_i, _ = self.blocks[b_id]
                self.bitmap_objects[vbo_i].data = new_bitmap[sl]

    def _load_textures(self):
        for block_id in sorted(self.blocks):
            vbo_i, block = self.blocks[block_id]
            n_data = np.abs(block.my_data[0]).copy(order="F").astype("float32").d
            n_data = (n_data - self.min_val) / (
                (self.max_val - self.min_val)
            )  # * self.diagonal)
            data_tex = Texture3D(data=n_data)
            bitmap_tex = Texture3D(
                data=block.source_mask * 255, min_filter="nearest", mag_filter="nearest"
            )
            self.texture_objects[vbo_i] = data_tex
            self.bitmap_objects[vbo_i] = bitmap_tex


class RGBAData(SceneData):
    name = "rgba_data"
    colormap_texture = traitlets.Instance(Texture1D)

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        va = VertexArray(name="tri", each=6)
        fq = FULLSCREEN_QUAD.reshape((6, 3), order="C")
        va.attributes.append(VertexAttribute(name="vertexPosition_modelspace", data=fq))
        return va

    def add_data(self, lines):
        assert lines.shape[1] == 4
        self.colormap_texture = Texture1D(boundary_x="clamp", data=lines)


class RGBADisplay(SceneComponent):
    name = "rgba_display"
    data = traitlets.Instance(RGBAData)
    priority = 20

    def _set_uniforms(self, scene, shader_program):
        shader_program._set_uniform("cm_texture", 0)

    def draw(self, scene, program):
        with self.data.colormap_texture.bind(0):
            GL.glDrawArrays(GL.GL_TRIANGLES, 0, 6)


class RGBALinePlot(SceneComponent):
    name = "rgba_line_plot"
    data = traitlets.Instance(LineData)
    priority = 20

    def draw(self, scene, program):
        for i, _channel in enumerate("rgba"):
            program._set_uniform("channel", i)
            GL.glDrawArrays(GL.GL_LINE_STRIP, 0, self.data.n_vertices)

    def _set_uniforms(self, scene, shader_program):
        shader_program._set_uniform(
            "viewport", np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype="f4")
        )


_cmaps = ["arbre", "viridis", "magma", "doom"]


class BlockRendering(SceneComponent):
    """
    A class that renders block data.  It may do this in one of several ways,
    including mesh outline.  This allows us to render a single collection of
    blocks multiple times in a single scene and to separate out the memory
    handling from the display.
    """

    name = "block_rendering"
    data = traitlets.Instance(BlockCollection)
    box_width = traitlets.CFloat(0.1)
    transfer_function = traitlets.Instance(TransferFunctionTexture)
    tf_min = traitlets.CFloat(0.0)
    tf_max = traitlets.CFloat(1.0)
    tf_log = traitlets.Bool(True)

    priority = 10

    def render_gui(self, imgui, renderer):
        changed = super(BlockRendering, self).render_gui(imgui, renderer)
        _, self.cmap_log = imgui.checkbox("Take log", self.cmap_log)
        changed = changed or _
        _, cmap_index = imgui.listbox(
            "Colormap", _cmaps.index(self.colormap.colormap_name), _cmaps
        )
        if _:
            self.colormap.colormap_name = _cmaps[cmap_index]
        changed = changed or _
        # Now, shaders
        shader_combos = list(sorted(component_shaders[self.name]))
        descriptions = [
            component_shaders[self.name][_]["description"] for _ in shader_combos
        ]
        selected = shader_combos.index(self.render_method)
        _, shader_ind = imgui.listbox("Shader", selected, descriptions)
        if _:
            self.render_method = shader_combos[shader_ind]
        changed = changed or _

        # Now for the transfer function stuff
        imgui.image_button(
            self.transfer_function.texture_name, 256, 32, frame_padding=0
        )
        imgui.text("Right click and drag to change")
        update = False
        data = self.transfer_function.data.astype("f4") / 255
        for i, c in enumerate("rgba"):
            imgui.plot_lines(
                f"## {c}",
                data[:, 0, i].copy(),
                scale_min=0.0,
                scale_max=1.0,
                graph_size=(256, 32),
            )
            if imgui.is_item_hovered() and imgui.is_mouse_dragging(2):
                update = True
                dx, dy = renderer.io.mouse_delta
                dy = -dy
                mi = imgui.get_item_rect_min()
                ma = imgui.get_item_rect_max()
                x, y = renderer.io.mouse_pos
                x = x - mi.x
                y = (ma.y - mi.y) - (y - mi.y)
                xb1 = floor(min(x + dx, x) * data.shape[0] / (ma.x - mi.x))
                xb2 = ceil(max(x + dx, x) * data.shape[0] / (ma.x - mi.x))
                yv1 = y / (ma.y - mi.y)
                yv2 = (y + dy) / (ma.y - mi.y)
                yv1, yv2 = (max(min(_, 1.0), 0.0) for _ in (yv1, yv2))
                if dx < 0:
                    yv2, yv1 = yv1, yv2
                    xb1 -= 1
                elif dx > 0:
                    xb2 += 1
                xb1 = max(0, xb1)
                xb2 = min(255, xb2)
                if renderer.io.key_shift:
                    yv1 = yv2 = 1.0
                elif renderer.io.key_ctrl:
                    yv1 = yv2 = 0.0
                data[xb1:xb2, 0, i] = np.mgrid[yv1 : yv2 : (xb2 - xb1) * 1j]
        if update:
            self.transfer_function.data = (data * 255).astype("u1")

    @traitlets.default("transfer_function")
    def _default_transfer_function(self):
        tf = TransferFunctionTexture(data=np.ones((256, 1, 4), dtype="u1") * 255)
        return tf

    def draw(self, scene, program):
        each = self.data.vertex_array.each
        with self.transfer_function.bind(target=2):
            for tex_ind, tex, bitmap_tex in self.data.viewpoint_iter(scene.camera):
                with tex.bind(target=0):
                    with bitmap_tex.bind(target=1):
                        GL.glDrawArrays(GL.GL_TRIANGLES, tex_ind * each, each)

    def _set_uniforms(self, scene, shader_program):
        cam = scene.camera
        shader_program._set_uniform("projection", cam.projection_matrix)
        shader_program._set_uniform("modelview", cam.view_matrix)
        shader_program._set_uniform(
            "viewport", np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype="f4")
        )
        shader_program._set_uniform("camera_pos", cam.position)
        shader_program._set_uniform("box_width", self.box_width)
        shader_program._set_uniform("ds_tex", 0)
        shader_program._set_uniform("bitmap_tex", 1)
        shader_program._set_uniform("tf_tex", 2)
        shader_program._set_uniform("tf_min", self.tf_min)
        shader_program._set_uniform("tf_max", self.tf_max)
        shader_program._set_uniform("tf_log", float(self.tf_log))


class BoxData(SceneData):
    name = "box_data"
    left_edge = YTPositionTrait([0.0, 0.0, 0.0])
    right_edge = YTPositionTrait([1.0, 1.0, 1.0])

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        va = VertexArray(name="box_outline", each=36)
        data = compute_box_geometry(self.left_edge, self.right_edge).copy()
        va.attributes.append(
            VertexAttribute(name="model_vertex", data=data.astype("f4"))
        )
        N = data.size // 4
        le = np.concatenate([[self.left_edge.copy()] for _ in range(N)])
        re = np.concatenate([[self.right_edge.copy()] for _ in range(N)])
        dds = self.right_edge - self.left_edge
        dds = np.concatenate([[dds.copy()] for _ in range(N)])
        va.attributes.append(VertexAttribute(name="in_left_edge", data=le.astype("f4")))
        va.attributes.append(
            VertexAttribute(name="in_right_edge", data=re.astype("f4"))
        )
        va.attributes.append(VertexAttribute(name="in_dx", data=dds.astype("f4")))
        return va


class BlockOutline(SceneAnnotation):
    """
    A class that renders outlines of block data.
    """

    name = "block_outline"
    data = traitlets.Instance(BlockCollection)
    box_width = traitlets.CFloat(0.1)
    box_color = traitlets.Tuple((1.0, 1.0, 1.0), trait=traitlets.CFloat())
    box_alpha = traitlets.CFloat(1.0)

    def draw(self, scene, program):
        each = self.data.vertex_array.each
        for tex_ind, _tex, _bitmap_tex in self.data.viewpoint_iter(scene.camera):
            GL.glDrawArrays(GL.GL_TRIANGLES, tex_ind * each, each)

    def render_gui(self, imgui, renderer):
        changed = super(BlockOutline, self).render_gui(imgui, renderer)
        _, bw = imgui.slider_float("Width", self.box_width, 0.001, 0.250)
        if _:
            self.box_width = bw
        changed = changed or _
        _, (r, g, b) = imgui.color_edit3("Color", *self.box_color)
        if _:
            self.box_color = (r, g, b)
        changed = changed or _
        _, ba = imgui.slider_float("Alpha", self.box_alpha, 0.0, 1.0)
        if _:
            self.box_alpha = ba
        changed = changed or _
        return changed

    def _set_uniforms(self, scene, shader_program):
        cam = scene.camera
        shader_program._set_uniform("projection", cam.projection_matrix)
        shader_program._set_uniform("modelview", cam.view_matrix)
        shader_program._set_uniform(
            "viewport", np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype="f4")
        )
        shader_program._set_uniform("camera_pos", cam.position)
        shader_program._set_uniform("box_width", self.box_width)
        shader_program._set_uniform("box_color", np.array(self.box_color))
        shader_program._set_uniform("box_alpha", self.box_alpha)


class BoxAnnotation(SceneAnnotation):
    name = "box_outline"
    data = traitlets.Instance(BoxData)
    box_width = traitlets.CFloat(0.05)
    box_alpha = traitlets.CFloat(1.0)

    def draw(self, scene, program):
        each = self.data.vertex_array.each
        GL.glDrawArrays(GL.GL_TRIANGLES, 0, each)

    def render_gui(self, imgui, renderer):
        changed = super(BoxAnnotation, self).render_gui(imgui, renderer)
        _, bw = imgui.slider_float("Width", self.box_width, 0.001, 0.250)
        if _:
            self.box_width = bw
        changed = changed or _
        _, ba = imgui.slider_float("Alpha", self.box_alpha, 0.0, 1.0)
        if _:
            self.box_alpha = ba
        changed = changed or _
        return changed

    def _set_uniforms(self, scene, shader_program):
        cam = scene.camera
        shader_program._set_uniform("projection", cam.projection_matrix)
        shader_program._set_uniform("modelview", cam.view_matrix)
        shader_program._set_uniform(
            "viewport", np.array(GL.glGetIntegerv(GL.GL_VIEWPORT), dtype="f4")
        )
        shader_program._set_uniform("camera_pos", cam.position)
        shader_program._set_uniform("box_width", self.box_width)
        shader_program._set_uniform("box_alpha", self.box_alpha)


class MeshData(SceneData):
    name = "mesh"
    data_source = traitlets.Instance(YTDataContainer)
    texture_objects = traitlets.Dict(trait=traitlets.Instance(Texture3D))
    texture_objects = traitlets.Dict(trait=traitlets.Instance(Texture3D))
    blocks = traitlets.Dict(default_value=())
    scale = traitlets.Bool(False)
    size = traitlets.CInt(-1)

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

        mesh = data_source.ds.index.meshes[mesh_id - 1]
        offset = mesh._index_offset
        vertices = mesh.connectivity_coords
        indices = mesh.connectivity_indices - offset

        data = data_source[field]

        return triangulate_mesh(vertices, data, indices)

    def add_data(self, field):
        v, d, i = self.get_mesh_data(self.data_source, field)
        v.shape = (v.size // 3, 3)
        d.shape = (d.size, 1)
        i.shape = (i.size, 1)
        i = i.astype("uint32")
        self.vertex_array.attributes.append(
            VertexAttribute(name="vertexPosition_modelspace", data=v)
        )
        self.vertex_array.attributes.append(VertexAttribute(name="vertexData", data=d))
        self.vertex_array.indices = i
        self.size = i.size

    @traitlets.default("vertex_array")
    def _default_vertex_array(self):
        return VertexArray(name="mesh_info", each=0)


class MeshRendering(SceneComponent):
    name = "mesh_rendering"
    data = traitlets.Instance(MeshData)

    def draw(self, scene, program):
        GL.glDrawElements(GL.GL_TRIANGLES, self.data.size, GL.GL_UNSIGNED_INT, None)

    def _set_uniforms(self, scene, shader_program):
        projection_matrix = scene.camera.projection_matrix
        view_matrix = scene.camera.view_matrix
        model_to_clip = np.dot(projection_matrix, view_matrix)
        shader_program._set_uniform("model_to_clip", model_to_clip)


class SceneGraph(traitlets.HasTraits):
    components = traitlets.List(
        trait=traitlets.Instance(SceneComponent), default_value=[]
    )
    annotations = traitlets.List(
        trait=traitlets.Instance(SceneAnnotation), default_value=[]
    )
    data_objects = traitlets.List(trait=traitlets.Instance(SceneData), default_value=[])
    camera = traitlets.Instance(IDVCamera)
    ds = traitlets.Instance(Dataset)
    fb = traitlets.Instance(Framebuffer, allow_none=True)
    input_captured_mouse = traitlets.Bool(False)
    input_captured_keyboard = traitlets.Bool(False)

    def add_volume(self, data_source, field_name, no_ghost=False):
        self.data_objects.append(BlockCollection(data_source=data_source))
        self.data_objects[-1].add_data(field_name, no_ghost=no_ghost)
        self.components.append(BlockRendering(data=self.data_objects[-1]))
        return self.components[-1]  # Only the rendering object

    def add_text(self, **kwargs):
        if "data" not in kwargs:
            if not any(_.name == "text_overlay" for _ in self.data_objects):
                self.data_objects.append(TextCharacters())
                self.data_objects[-1].build_textures()
            kwargs["data"] = next(
                (_ for _ in self.data_objects if _.name == "text_overlay")
            )
        self.annotations.append(TextAnnotation(**kwargs))
        return self.annotations[-1]

    def add_box(self, left_edge, right_edge):
        data = BoxData(left_edge=left_edge, right_edge=right_edge)
        self.data_objects.append(data)
        self.annotations.append(BoxAnnotation(data=data))

    def __iter__(self):
        elements = self.components + self.annotations
        for element in sorted(elements, key=lambda a: a.priority):
            yield element

    def render(self):
        origin_x, origin_y, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        with self.bind_buffer():
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for element in self:
            # If we need to clear the region, we need to use the scissor test
            db = element.display_bounds
            new_origin_x = origin_x + width * db[0]
            new_origin_y = origin_y + height * db[2]
            new_width = (db[1] - db[0]) * width
            new_height = (db[3] - db[2]) * height
            GL.glViewport(
                int(new_origin_x), int(new_origin_y), int(new_width), int(new_height)
            )
            if element.clear_region:
                GL.glEnable(GL.GL_SCISSOR_TEST)
                GL.glScissor(
                    int(new_origin_x),
                    int(new_origin_y),
                    int(new_width),
                    int(new_height),
                )
                GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
                GL.glDisable(GL.GL_SCISSOR_TEST)
            element.run_program(self)
        GL.glViewport(origin_x, origin_y, width, height)

    @contextlib.contextmanager
    def bind_buffer(self):
        if self.fb is not None:
            with self.fb.bind(clear=False):
                yield
        else:
            yield

    @property
    def image(self):
        if self.fb is not None:
            arr = self.fb.data[::-1, :, :]
            arr.swapaxes(0, 1)
            return arr
        _, _, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        arr = GL.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_FLOAT)[::-1, :, :]
        arr.swapaxes(0, 1)
        return arr

    @property
    def depth(self):
        if self.fb is not None:
            arr = self.fb.depth_data[::-1, :]
            arr.swapaxes(0, 1)
            return arr
        _, _, width, height = GL.glGetIntegerv(GL.GL_VIEWPORT)
        arr = GL.glReadPixels(0, 0, width, height, GL.GL_DEPTH_COMPONENT, GL.GL_FLOAT)[
            ::-1, :
        ]
        arr.swapaxes(0, 1)
        return arr

    @staticmethod
    def from_ds(ds, field, no_ghost=True):
        # Here we make a bunch of guesses.  Nothing too complex -- only two
        # arguments: dataset and field.  If you supply a rendering context,
        # great.  If not, we'll make one.
        if isinstance(ds, Dataset):
            data_source = ds.all_data()
        else:
            ds, data_source = ds.ds, ds
        center = ds.domain_center
        pos = center + 1.5 * ds.domain_width
        near_plane = 3.0 * ds.index.get_smallest_dx().min().d

        c = TrackballCamera(position=pos, focus=center, near_plane=near_plane)

        scene = SceneGraph(camera=c)
        scene.add_volume(data_source, field, no_ghost=no_ghost)
        return scene
