# encoding: utf-8
"""
Shader and ShaderProgram wrapper classes for vertex and fragment shaders used 
in Interactive Data Visualization
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
from contextlib import contextmanager
import traitlets
import numpy as np

try:
    from contextlib import ExitStack
except ImportError:
    try:
        from contextlib2 import ExitStack
    except ImportError:
        raise RuntimeError("IDV requires either Python 3.3+ or contextlib2")


# Set up a mapping from numbers to names

const_types = (GL.constant.IntConstant,
               GL.constant.LongConstant,
               GL.constant.FloatConstant)

num_to_const = {}
for i in dir(GL):
    if i.startswith("GL_"):
        v = getattr(GL, i)
        if not isinstance(v, const_types): continue
        num_to_const[v.real] = v

_coersion_funcs = {
        'FLOAT': float,
        'DOUBLE': float,
        'INT': int,
        'UNSIGNED': int,
        'BOOL': bool
}

_shapes = {
        'MAT2'   : (2, 2),
        'MAT3'   : (3, 3),
        'MAT4'   : (4, 4),
        'MAT2x3' : (2, 3),
        'MAT2x4' : (2, 4),
        'MAT3x2' : (3, 2),
        'MAT3x4' : (3, 4),
        'MAT4x2' : (4, 2),
        'MAT4x3' : (4, 3),
        'VEC2'   : (2,),
        'VEC3'   : (3,),
        'VEC4'   : (4,),
}

# PyOpenGL has the reverse mapping for this
gl_to_np = {
        'FLOAT'    : 'f',
        'DOUBLE'   : 'd',
        'INT'      : 'i',
        'UNSIGNED' : 'I',
        'BOOL'     : 'b',
}

def coerce_uniform_type(val, gl_type):
    # gl_type here must be in const_types
    if not isinstance(gl_type, const_types):
        gl_type = num_to_const[gl_type]
    # Now we can get down to business!
    spec = gl_type.name.split("_")[1:] # Strip out the GL_
    # We know what to do with:
    #   FLOAT DOUBLE INT UNSIGNED BOOL
    # We can ignore:
    #    SAMPLER IMAGE 
    if "SAMPLER" in spec or "IMAGE" in spec:
        # Do nothing to these, and let PyOpenGL handle it
        return val
    if len(spec) == 1 or spec == ['UNSIGNED', 'INT']:
        return _coersion_funcs[spec[0]](val)
    # We need to figure out if it's a matrix, a vector, etc.
    shape = _shapes[spec[-1]]
    dtype = gl_to_np[spec[0]]
    val = np.asanyarray(val, dtype = dtype)
    val.shape = shape
    return val

class TextureBoundary(traitlets.TraitType):
    default_value = GL.GL_CLAMP_TO_EDGE
    info_text = "A boundary type of mirror, clamp, or repeat"

    def validate(self, obj, value):
        if isinstance(value, str):
            try:
                return {'clamp': GL.GL_CLAMP_TO_EDGE,
                        'mirror': GL.GL_MIRRORED_REPEAT,
                        'repeat': GL.GL_REPEAT}[value.lower()]
            except KeyError:
                self.error(obj, value)
        elif value in (GL.GL_CLAMP_TO_EDGE,
                       GL.GL_MIRRORED_REPEAT,
                       GL.GL_REPEAT):
            return value
        self.error(obj, value)

TEX_CHANNELS = {
        1: (GL.GL_R32F, GL.GL_RED),
        2: (GL.GL_RG32F, GL.GL_RG),
        3: (GL.GL_RGB32F, GL.GL_RGB),
        4: (GL.GL_RGBA32F, GL.GL_RGBA)
}

class Texture(traitlets.HasTraits):
    texture_name = traitlets.CInt(-1)
    data = traitlets.Instance(np.ndarray)

    @traitlets.default('texture_name')
    def _default_texture_name(self):
        return GL.glGenTextures(1)

    @contextmanager
    def bind(self):
        GL.glBindTexture(self.dim_enum, self.texture_name)
        yield
        GL.glBindTexture(self.dim_enum, 0)

class Texture1D(Texture):
    boundary_x = TextureBoundary()
    dims = 1
    dim_enum = GL.GL_TEXTURE_1D

    @traitlets.observe("data")
    def _set_data(self, change):
        with self.bind():
            data = change['new']
            if len(data.shape) == 3:
                channels = data.shape[-1]
            else:
                channels = 1
            dx = data.shape[0]
            type1, type2 = TEX_CHANNELS[channels]
            GL.glTexStorage1D(GL.GL_TEXTURE_1D, 1, type1, dx)
            GL.glTexSubImage1D(GL.GL_TEXTURE_1D, 0, 0, dx,
                        type2, GL.GL_FLOAT, data)
            GL.glTexParameterf(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_WRAP_S,
                    self.boundary_x)
            GL.glGenerateMipmap(GL.GL_TEXTURE_1D)

class Texture2D(Texture):
    boundary_x = TextureBoundary()
    boundary_y = TextureBoundary()
    dims = 2
    channels = 1
    dim_enum = GL.GL_TEXTURE_2D

    @traitlets.observe("data")
    def _set_data(self, change):
        with self.bind():
            data = change['new']
            if len(data.shape) == 3:
                channels = data.shape[-1]
            else:
                channels = 1
            dx, dy = data.shape[:2]
            type1, type2 = TEX_CHANNELS[channels]
            GL.glTexStorage2D(GL.GL_TEXTURE_2D, 1, type1, dx, dy)
            GL.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, dx, dy, 
                        type2, GL.GL_FLOAT, data)
            GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S,
                    self.boundary_x)
            GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T,
                    self.boundary_y)
            GL.glGenerateMipmap(GL.GL_TEXTURE_2D)

class Texture3D(Texture):
    boundary_x = TextureBoundary()
    boundary_y = TextureBoundary()
    boundary_z = TextureBoundary()
    dims = 3
    dim_enum = GL.GL_TEXTURE_3D

    @traitlets.observe("data")
    def _set_data(self, change):
        with self.bind():
            data = change['new']
            if len(data.shape) == 3:
                channels = data.shape[-1]
            else:
                channels = 1
            dx, dy, dz = data.shape[:3]
            type1, type2 = TEX_CHANNELS[channels]
            GL.glTexStorage3D(GL.GL_TEXTURE_3D, 1, type1, dx, dy, dz)
            GL.glTexSubImage3D(GL.GL_TEXTURE_3D, 0, 0, 0, 0, dx, dy, dz,
                        type2, GL.GL_FLOAT, data)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S,
                    self.boundary_x)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T,
                    self.boundary_y)
            GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R,
                    self.boundary_z)
            GL.glGenerateMipmap(GL.GL_TEXTURE_3D)

class VertexAttribute(traitlets.HasTraits):
    name = traitlets.CUnicode("attr")
    id = traitlets.CInt(-1)
    data = traitlets.Instance(np.ndarray)
    each = traitlets.CInt(-1)

    @traitlets.default('id')
    def _id_default(self):
        return GL.glGenBuffers(1)

    @contextmanager
    def bind(self, program = None):
        loc = -1
        if program is not None:
            loc = GL.glGetAttribLocation(program, self.name)
            if loc < 0:
                return -1
            GL.glEnableVertexAttribArray(loc)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.id)
        if loc >= 0:
            GL.glVertexAttribPointer(loc, self.each, GL.GL_FLOAT, False, 0,
                    None)
        yield
        if loc >= 0:
            GL.glDisableVertexAttribArray(loc)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

    @traitlets.observe("data")
    def _set_data(self, change):
        arr = change['new']
        self.each = arr.shape[-1]
        with self.bind():
            GL.glBufferData(GL.GL_ARRAY_BUFFER, arr.nbytes, arr, GL.GL_STATIC_DRAW)

class VertexArray(traitlets.HasTraits):
    name = traitlets.CUnicode("vertex")
    id = traitlets.CInt(-1)
    element_ids = traitlets.Instance(np.ndarray, allow_none = True)
    attributes = traitlets.List(trait=traitlets.Instance(VertexAttribute))
    each = traitlets.CInt(-1)

    @traitlets.default('id')
    def _id_default(self):
        return GL.glGenVertexArrays(1)

    @contextmanager
    def bind(self, program = None):
        GL.glBindVertexArray(self.id)
        # We only bind the attributes if we have a program too
        if program is None:
            attrs = []
        else:
            attrs = self.attributes
        with ExitStack() as stack:
            _ = [stack.enter_context(_) for _ in attrs]
            yield
        GL.glBindVertexArray(0)

    @traitlets.observe("element_ids")
    def _set_elements(self, change):
        arr = change['new']

class Framebuffer(traitlets.HasTraits):
    fb_id = traitlets.CInt(-1)
    db_id = traitlets.CInt(-1)
    fb_tex = traitlets.Instance(Texture2D)
    viewport = traitlets.Tuple(
            traitlets.Int(), traitlets.Int(),
            traitlets.Int(), traitlets.Int())
    initialized = traitlets.Bool(False)

    @property
    def data(self):
        _, _, width, height = self.viewport
        with self.bind(clear = False):
            debug_buffer = GL.glReadPixels(0, 0, width, height, GL.GL_RGB,
                                           GL.GL_UNSIGNED_BYTE)
            arr = np.fromstring(debug_buffer, "uint8", count = width*height*3)
        return arr.reshape((width, height, 3))

    @traitlets.default("viewport")
    def _viewport_default(self):
        # origin_x, origin_y, width, height
        return GL.glGetIntegerv(GL.GL_VIEWPORT)

    @traitlets.default("fb_id")
    def _fb_id_default(self):
        return GL.GenFrameBuffers(1)

    @traitlets.default("db_id")
    def _db_id_default(self):
        return GL.GenRenderBuffers(1)

    def _fb_tex_default(self):
        data = np.zeros( (self.viewport[2], self.viewport[3], 4), "f4")
        return Texture2D(data = data, boundary_x = "repeat", 
                boundary_y = "repeat")

    @contextmanager
    def bind(self, clear = True):
        if not self.initialized:
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fb_id)
            GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, self.db_id)
            GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT32F,
                                     self.viewport[2], self.viewport[3])
            GL.glFramebufferRenderbuffer(
                GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT, GL.GL_RENDERBUFFER,
                self.db_id
            )

            GL.glFramebufferTexture2D(
                GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D,
                self.fb_tex,
                0 # mipmap level, normally 0
            )
            status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
            if status != GL.GL_FRAMEBUFFER_COMPLETE:
                raise RuntimeError
        if clear:
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fb_id)
        yield
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
