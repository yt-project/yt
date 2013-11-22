"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import OpenGL.GL as GL
import OpenGL.GLUT as GLUT
import OpenGL.GLU as GLU
import OpenGL.GL.shaders as shaders
from OpenGL.arrays import vbo, ArrayDatatype
import OpenGL.GL.ARB.framebuffer_object as GL_fbo
from functools import wraps

def render_fbo(func):
    @wraps(func)
    def fbo_renderer(self, fbo_id = "fbo_id", fbo_depth = "fbo_depth", **kwargs):
        GL_fbo.glBindFramebuffer(GL_fbo.GL_FRAMEBUFFER,
                                 self.gl_state[fbo_id])
        GL_fbo.glBindRenderbuffer(GL_fbo.GL_RENDERBUFFER,
                                  self.gl_state[fbo_depth])
        self.set_viewport()
        GL.glPushAttrib(GL.GL_VIEWPORT_BIT)
        GL.glViewport(0, 0, self.gl_state["width"], self.gl_state["height"])
        status = GL_fbo.glCheckFramebufferStatus(GL_fbo.GL_FRAMEBUFFER)
        assert(status == GL_fbo.GL_FRAMEBUFFER_COMPLETE)

        func(self, **kwargs)

        GL.glPopAttrib(GL.GL_VIEWPORT_BIT)
        GL_fbo.glBindFramebuffer(GL_fbo.GL_FRAMEBUFFER, 0)

    return fbo_renderer

def create_fbo(gl_state, fbo_id = "fbo_id", fbo_depth = "fbo_depth",
               fbo_tex = "fbo_tex"):
    GL.glActiveTexture(GL.GL_TEXTURE0)
    id_fbo = GL.glGenTextures(1)
    id_depth = GL_fbo.glGenRenderbuffers(1)
    gl_state[fbo_tex] = id_fbo
    gl_state[fbo_depth] = id_depth
    gl_state[fbo_id] = GL_fbo.glGenFramebuffers(1)

    GL.glBindTexture(GL.GL_TEXTURE_2D, id_fbo)
    GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8,
                    gl_state["width"], gl_state["height"], 0,
                    GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, None)

    GL_fbo.glBindFramebuffer(GL_fbo.GL_FRAMEBUFFER, gl_state[fbo_id])
    GL_fbo.glFramebufferTexture2D(GL_fbo.GL_FRAMEBUFFER,
        GL_fbo.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D, gl_state[fbo_tex], 0)

    GL_fbo.glBindRenderbuffer(GL_fbo.GL_RENDERBUFFER, id_depth)
    GL_fbo.glRenderbufferStorage(GL_fbo.GL_RENDERBUFFER,
        GL.GL_DEPTH_COMPONENT, gl_state["width"], gl_state["height"])
    GL_fbo.glFramebufferRenderbuffer(
        GL_fbo.GL_FRAMEBUFFER, GL_fbo.GL_DEPTH_ATTACHMENT,
        GL_fbo.GL_RENDERBUFFER, id_depth)

    GL_fbo.glBindFramebuffer(GL_fbo.GL_FRAMEBUFFER, 0)
    GL_fbo.glBindRenderbuffer(GL_fbo.GL_RENDERBUFFER, 0)

def identity_view(func):
    pass

def translate_view(func):
    pass
