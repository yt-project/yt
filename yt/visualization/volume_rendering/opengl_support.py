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
