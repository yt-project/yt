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

import os
import OpenGL.GL as GL
import contextlib
from yt.extern.six import add_metaclass
from collections import OrderedDict
from yt.utilities.exceptions import \
    YTInvalidShaderType, \
    YTUnknownUniformKind, \
    YTUnknownUniformSize

known_shaders = {}

class ShaderProgram(object):
    '''
    Wrapper class that compiles and links vertex and fragment shaders
    into a shader program.

    Parameters
    ----------

    vertex_shader : string or :class:`yt.visualization.volume_rendering.shader_objects.VertexShader`
        The vertex shader used in the Interactive Data Visualization pipeline.

    fragment_shader : string or :class:`yt.visualization.volume_rendering.shader_objects.FragmentShader`
        The fragment shader used in the Interactive Data Visualization pipeline.
    '''

    def __init__(self, vertex_shader=None, fragment_shader=None):
        # Don't allow just one.  Either neither or both.
        if vertex_shader is None and fragment_shader is None:
            pass
        elif None not in (vertex_shader, fragment_shader):
            self.link(vertex_shader, fragment_shader)
        else:
            raise RuntimeError
        self._uniform_funcs = OrderedDict()

    def link(self, vertex_shader, fragment_shader):
        # There are more types of shaders, but for now we only allow v&f.
        self.program = GL.glCreateProgram()
        if not isinstance(vertex_shader, VertexShader):
            vertex_shader = VertexShader(vertex_shader)
        if not isinstance(fragment_shader, FragmentShader):
            fragment_shader = FragmentShader(fragment_shader)
        self.vertex_shader = vertex_shader
        self.fragment_shader = fragment_shader
        GL.glAttachShader(self.program, vertex_shader.shader)
        GL.glAttachShader(self.program, fragment_shader.shader)
        GL.glLinkProgram(self.program)
        result = GL.glGetProgramiv(self.program, GL.GL_LINK_STATUS)
        if not result:
            raise RuntimeError(GL.glGetProgramInfoLog(self.program))
        vertex_shader.delete_shader()
        fragment_shader.delete_shader()

    def delete_program(self):
        if self.program is not None:
            GL.glDeleteProgram(self.program)
            self.program = None

    def _guess_uniform_func(self, value):
        # We make a best-effort guess.
        # This does NOT work with arrays of uniforms.
        # First look at the dtype kind.  Fortunately, this falls into either
        # 'f' or 'i', which matches nicely with OpenGL.
        # Note that in some implementations, it seems there is also a 'd' type,
        # but we will not be using that here.
        if isinstance(value, int):
            return GL.glUniform1i
        elif isinstance(value, float):
            return GL.glUniform1f
        else:
            kind = value.dtype.kind
        if kind not in 'if':
            raise YTUnknownUniformKind(kind)
        if len(value.shape) == 1:
            if value.size > 4:
                raise YTUnknownUniformSize(value.size)
            func = self._set_scalar_uniform(kind, value.size)
        elif len(value.shape) == 2:
            if value.shape[0] != value.shape[1]:
                raise YTUnknownUniformSize(value.shape)
            func = self._set_matrix_uniform(kind, value.shape)
        else:
            raise YTUnknownUniformSize(value.shape)
        return func

    def _set_scalar_uniform(self, kind, size_spec):
        gl_func = getattr(GL, "glUniform%s%sv" % (size_spec, kind))
        def _func(location, value):
            return gl_func(location, 1, value)
        return _func

    def _set_matrix_uniform(self, kind, size_spec):
        assert(size_spec[0] == size_spec[1])
        gl_func = getattr(GL, "glUniformMatrix%s%sv" % (size_spec[0], kind))
        def _func(location, value):
            return gl_func(location, 1, GL.GL_TRUE, value)
        return _func

    def _set_uniform(self, name, value):
        # We need to figure out how to pass it in.
        if name not in self._uniform_funcs:
            self._uniform_funcs[name] = self._guess_uniform_func(value)
        loc = GL.glGetUniformLocation(self.program, name)
        return self._uniform_funcs[name](loc, value)

    @contextlib.contextmanager
    def enable(self):
        GL.glUseProgram(self.program)
        yield
        GL.glUseProgram(0)

    def bind_vert_attrib(self, name, bind_loc, size):
        loc = GL.glGetAttribLocation(self.program, name)
        GL.glEnableVertexAttribArray(loc)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, bind_loc)
        GL.glVertexAttribPointer(loc, size, GL.GL_FLOAT, False, 0, None)

    def disable_vert_attrib(self, name):
        loc = GL.glGetAttribLocation(self.program, name)
        GL.glDisableVertexAttribArray(loc)

class RegisteredShader(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if getattr(cls, "_shader_name", None) is not None:
            known_shaders[cls._shader_name] = cls

@add_metaclass(RegisteredShader)
class Shader(object):
    '''
    Creates a shader from source

    Parameters
    ----------

    source : str
        This can either be a string containing a full source of a shader,
        an absolute path to a source file or a filename of a shader
        residing in the ./shaders/ directory.

    '''

    _shader = None
    _source = None
    _shader_name = None

    def __init__(self, source=None):
        if source:
            self.compile(source)

    def _get_source(self, source):
        if ";" in source:
            # This is probably safe, right?  Enh, probably.
            return source
        if os.path.isfile(source):
            sh_directory = ''
        else:
            sh_directory = os.path.join(os.path.dirname(__file__), "shaders")
        fn = os.path.join(sh_directory, source)
        if not os.path.isfile(fn):
            raise YTInvalidShaderType(source)
        return open(fn, 'r').read()

    def compile(self, source = None, parameters = None):
        if source is None:
            source = self._source
            if source is None: raise RuntimeError
        if parameters is not None:
            raise NotImplementedError
        source = self._get_source(source)
        shader_type_enum = getattr(GL,
            'GL_%s_SHADER' % self.shader_type.upper())
        shader = GL.glCreateShader(shader_type_enum)
        # We could do templating here if we wanted.
        self.shader_source = source
        GL.glShaderSource(shader, source)
        GL.glCompileShader(shader)
        result = GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS)
        if not(result):
            raise RuntimeError(GL.glGetShaderInfoLog(shader))
        self._shader = shader

    @property
    def shader(self):
        if self._shader is None:
            self.compile()
        return self._shader

    def delete_shader(self):
        if self.shader is not None:
            GL.glDeleteShader(self.shader)
            self._shader = None

    def __del__(self):
        # This is not guaranteed to be called
        self.delete_shader()

class FragmentShader(Shader):
    '''Wrapper class for fragment shaders'''
    shader_type = "fragment"

class VertexShader(Shader):
    '''Wrapper class for vertex shaders'''
    shader_type = "vertex"

class ApplyColormapFragmentShader(FragmentShader):
    '''A second pass fragment shader used to apply a colormap to the result of
    the first pass rendering'''
    _source = "apply_colormap.fragmentshader"
    _shader_name = "apply_colormap.f"

class MaxIntensityFragmentShader(FragmentShader):
    '''A first pass fragment shader that computes Maximum Intensity Projection
    of the data. See :ref:`projection-types` for more information.'''
    _source = "max_intensity.fragmentshader"
    _shader_name = "max_intensity.f"

class NoOpFragmentShader(FragmentShader):
    '''A second pass fragment shader that performs no operation. Usually used if
    the first pass already took care of applying proper color to the data'''
    _source = "noop.fragmentshader"
    _shader_name = "noop.f"

class PassthroughFragmentShader(FragmentShader):
    '''A first pass fragment shader that performs no operation. Used for debug
    puproses. It's distinct from NoOpFragmentShader, because of the number of
    uniforms'''
    _source = "passthrough.fragmentshader"
    _shader_name = "passthrough.f"

class ProjectionFragmentShader(FragmentShader):
    '''A first pass fragment shader that performs unweighted integration of the
    data along the line of sight. See :ref:`projection-types` for more
    information.'''
    _source = "projection.fragmentshader"
    _shader_name = "projection.f"

class TransferFunctionFragmentShader(FragmentShader):
    '''A first pass fragment shader that performs ray casting using transfer
    function. See :ref:`volume-rendering-method` for more details.'''
    _source = "transfer_function.fragmentshader"
    _shader_name = "transfer_function.f"

class DefaultVertexShader(VertexShader):
    '''A first pass vertex shader that tranlates the location of vertices from
    the world coordinates to the viewing plane coordinates'''
    _source = "default.vertexshader"
    _shader_name = "default.v"

class PassthroughVertexShader(VertexShader):
    '''A second pass vertex shader that performs no operations on vertices'''
    _source = "passthrough.vertexshader"
    _shader_name = "passthrough.v"

class MeshVertexShader(VertexShader):
    '''A vertex shader used for unstructured mesh rendering.'''
    _source = "mesh.vertexshader"
    _shader_name = "mesh.v"

class MeshFragmentShader(FragmentShader):
    '''A vertex shader used for unstructured mesh rendering.'''
    _source = "mesh.fragmentshader"
    _shader_name = "mesh.f"
