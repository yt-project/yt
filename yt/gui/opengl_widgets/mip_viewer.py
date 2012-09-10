"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation:  UCSD
License:
  Copyright (C) 2010-2011 Matthew Turk  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# A simple image viewer, soon to be useful for stereo images, using OpenGL.

import sys, os
import OpenGL.GL as GL
import OpenGL.GLUT as GLUT
import OpenGL.GLU as GLU
import OpenGL.GL.shaders as shaders
from OpenGL.arrays import vbo, ArrayDatatype
import OpenGL.GL.ARB.framebuffer_object as GL_fbo
import Image
import glob
import numpy as np
import time

from small_apps import ViewHandler3D, GenericGLUTScene
from yt.visualization.image_writer import map_to_colors

from rendering_contexts import render_fbo, create_fbo, identity_view, \
        translate_view
from yt.visualization.volume_rendering.api import \
    HomogenizedVolume

ESCAPE = '\033'

_verts = ( (0,0,0), (1,0,0), (1,1,0), (0,1,0),
           (0,0,0), (1,0,0), (1,0,1), (0,0,1),
           (0,0,0), (0,0,1), (0,1,1), (0,1,0),
           (0,1,0), (1,1,0), (1,1,1), (0,1,1),
           (1,1,0), (1,0,0), (1,0,1), (1,1,1),
           (0,0,1), (0,1,1), (1,1,1), (1,0,1) )

_verts = ( (1,1,0), (0,1,0), (0,1,1), (1,1,1),
           (1,0,1), (0,0,1), (0,0,0), (1,0,0),
           (1,1,1), (0,1,1), (0,0,1), (1,0,1),
           (1,0,0), (0,0,0), (0,1,0), (1,1,0),
           (0,1,1), (0,1,0), (0,0,0), (0,0,1),
           (1,1,0), (1,1,1), (1,0,1), (1,0,0) )

_corner_list = [0,1,2,3, 4,5,6,7, 3,2,6,5, 0,4,7,1, 0,3,5,4, 1,7,6,2]
_corner_vals = [ (0,0,0), (1,0,0), (1,1,0), (0,1,0),
                 (0,0,1), (0,1,1), (1,1,1), (1,0,1) ]

_reversed = {2:0,1:1,0:2}

def _compress(lr, s, i):
    
    if lr == 0: return 0.5/s[_reversed[i]]
    else: return -0.5/s[_reversed[i]]
    

class MIPScene(GenericGLUTScene):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH)
    _title = "MIP"

    gl_state = None

    def _get_brick_vertices(self, offset):
        for b in self.hv.bricks:
            s = [ [-0.5, -0.5, -0.5],
                  [b.my_data[0].shape[i] - 1.5 for i in reversed(xrange(3))] ]
            for corner in _corner_list:
                for i, v in enumerate(_corner_vals[corner]):
                    yield s[v][i]

    def _get_texture_vertices(self):
        vs = [np.zeros(3, dtype='float32'),
              np.ones(3, dtype='float32')]
        #vs.reverse()
        for b in self.hv.bricks:
            shape = b.my_data[0].shape
            for corner in _corner_list:
                for i,v in enumerate(_corner_vals[corner]):
                    yield vs[v][i] + _compress(v, shape, i)

    def _setup_bricks(self):
        self._brick_textures = []
        for g in self.hv.bricks:
            self._upload_brick_textures(g)

    def _upload_brick_textures(self, brick):
        ix, iy, iz = brick.my_data[0].shape

        GL.glActiveTexture(GL.GL_TEXTURE0)
        id_field = GL.glGenTextures(1)
        upload = brick.my_data[0].astype("float32")
        #upload = (upload - -31.847) / ( -25.948 - -31.847 )
        #mi, ma = -31.847, -25.948
        #mi, ma = -27.2062, -20.9649 
        #upload = (upload - mi)/(ma - mi)
        self.mi = min(upload.min(), self.mi)
        self.ma = max(upload.max(), self.ma)
        #upload = (255*(upload - -31.0) / (-25.0 - -31.0)).astype("uint8")
        
        GL.glActiveTexture(GL.GL_TEXTURE0)
        GL.glBindTexture(GL.GL_TEXTURE_3D, id_field)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexImage3D(GL.GL_TEXTURE_3D, 0, GL.GL_LUMINANCE32F_ARB, iz, iy, ix, 0,
                        GL.GL_LUMINANCE, GL.GL_FLOAT, upload)

        DW = self.hv.pf.domain_right_edge - self.hv.pf.domain_left_edge
        dds = ((brick.RightEdge - brick.LeftEdge) /
               (np.array([ix,iy,iz], dtype='float32')-1)) / DW
        BLE = brick.LeftEdge / DW - 0.5
        self._brick_textures.append(
            (id_field, (ix-1,iy-1,iz-1), dds, BLE))

        print "Uploaded", len(self._brick_textures)

    def _setup_colormap(self):

        buffer = np.mgrid[0.0:1.0:256j]
        colors = map_to_colors(buffer, "algae")
        
        GL.glActiveTexture(GL.GL_TEXTURE1)
        id_cmap = GL.glGenTextures(1)

        GL.glBindTexture(GL.GL_TEXTURE_1D, id_cmap)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexParameterf(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameterf(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameterf(GL.GL_TEXTURE_1D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexImage1D(GL.GL_TEXTURE_1D, 0, GL.GL_RGBA, 256, 0,
                        GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, colors)

        self.gl_state["cmap_tex"] = id_cmap

    def __init__(self, hv, offset = 0.5):
        self.offset = offset
        self.mi, self.ma = 1e30, -1e30
        self.hv = hv
        self.coord = 0.0
        self.tfac = 10.0
        self.rfac = 0.5
        self.wireframe = False
        self.glda = True
        self._setup_keypress_handler()
        self.gl_state = {}
        GenericGLUTScene.__init__(self, 800, 800)

        num = len(hv.bricks) * 6 * 4
        self.v = np.fromiter(self._get_brick_vertices(offset),
                             dtype = 'float32', count = num * 3)
        self.vertices = vbo.VBO(self.v)

        self.t = np.fromiter(self._get_texture_vertices(),
                             dtype = 'float32', count = num * 3)
        self.tvertices = vbo.VBO(self.t)

        self.ng = len(hv.bricks)
        self.position = np.zeros(3, dtype='float')
        self.rotation = np.zeros(3, dtype='float') + 30
        self.position[2] = -2 # Offset backwards a bit

        self._setup_bricks()
        create_fbo(self.gl_state)
        self._setup_colormap()

    def init_opengl(self, width, height):
        # One-time GL setup

        self.gl_state["width"] = width
        self.gl_state["height"] = height

        self.set_viewport()

        # Now we compile our shaders
        self.recompile()

    def set_viewport(self):
        if self.wireframe:
            GL.glClearColor(1, 1, 1, 1)
        else:
            GL.glClearColor(0, 0, 0, 1)
        GL.glColor3f(1, 0, 0)
        GL.glEnable(GL.GL_DEPTH_TEST)

        GL.glViewport(0, 0, self.gl_state["width"], self.gl_state["height"])
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GLU.gluPerspective(60., self.gl_state["width"] / float(self.gl_state["height"]), 1e-3, 10.)
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def recompile(self):
        base = os.path.dirname(__file__) + "/"
        self.program = compileProgram(
            shaders.compileShader(
                open(base+"calculateRay.vertex.glsl").read(),
                GL.GL_VERTEX_SHADER),
            shaders.compileShader(
                open(base+"mip.fragment.glsl").read(),
                GL.GL_FRAGMENT_SHADER))

        def _set_simple_uniform(prog):
            GL.glUseProgram(prog)
            loc = GL.glGetUniformLocation(prog, 'buffer')
            GL.glUniform1i(loc, 0)
            loc = GL.glGetUniformLocation(prog, 'colormap')
            GL.glUniform1i(loc, 1)
            GL.glUseProgram(0)

        self.fbo_program = compileProgram(
            shaders.compileShader(
                open(base+"framebuffer.vertex.glsl").read(),
                GL.GL_VERTEX_SHADER),
            shaders.compileShader(
                open(base+"colormap.fragment.glsl").read(),
                GL.GL_FRAGMENT_SHADER),
            callback = _set_simple_uniform)

        self.uniform_locations = dict( (
                (v, GL.glGetUniformLocation(self.program, v)) for v in
                ['texture','shape','scaleBias','stepRatio','position']
            ) )
        self.fbo_uniform_locations = dict( (
                (v, GL.glGetUniformLocation(self.fbo_program, v)) for v in
                ['colormap','buffer']
            ) )

        print self.uniform_locations, self.fbo_uniform_locations

    @render_fbo
    def _draw_boxes(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        if self.wireframe:
            GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE)
        else:
            GL.glBlendEquation(GL.GL_MAX)
            GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)
            GL.glDisable(GL.GL_LINE_STIPPLE)
            GL.glDisable(GL.GL_LINE_SMOOTH)
            GL.glDisable(GL.GL_POINT_SMOOTH)
            GL.glEnable(GL.GL_CULL_FACE)
            #GL.glCullFace(GL.GL_FRONT)
            GL.glCullFace(GL.GL_BACK)

        GL.glLoadIdentity()
        GL.glTranslatef(*self.position)
        GL.glRotatef(self.rotation[0], 0, 0, 1)
        GL.glRotatef(self.rotation[1], 0, 1, 0)
        GL.glRotatef(self.rotation[2], 1, 0, 0)

        GL.glColor3f(0.0, 0.0, 0.0)
        GL.glUseProgram(self.program)
        scalebias = ( 1.0/(self.ma - self.mi), -self.mi)
        GL.glUniform1i(self.uniform_locations["texture"], 0)
        GL.glUniform2f(self.uniform_locations["scaleBias"], scalebias[0], scalebias[1])
        GL.glUniform1f(self.uniform_locations["stepRatio"], 1.0)
        GL.glUniform3f(self.uniform_locations["position"], *self.position)

        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glEnableClientState(GL.GL_TEXTURE_COORD_ARRAY)
        #GL.glEnableClientState(GL.GL_TEXTURE_COORD_ARRAY_EXT)

        GL.glColor3f(0.0, 0.0, 0.0)

        GL.glActiveTexture(GL.GL_TEXTURE0)
        GL.glEnable(GL.GL_TEXTURE_3D)

        self.tvertices.bind()
        GL.glTexCoordPointer(3, GL.GL_FLOAT, 0, self.tvertices)
        self.vertices.bind()
        GL.glVertexPointer(3, GL.GL_FLOAT, 0, self.vertices)

        for i, texinfo in enumerate(self._brick_textures):
            tex, shape, width, LE = texinfo
            GL.glActiveTexture(GL.GL_TEXTURE0)
            GL.glBindTexture(GL.GL_TEXTURE_3D, tex)
            GL.glUniform1f(self.uniform_locations["stepRatio"], 0.1)
            GL.glUniform3f(self.uniform_locations["shape"],
                shape[2], shape[1], shape[0])
            GL.glPushMatrix()
            GL.glTranslate(LE[2], LE[1], LE[0])
            GL.glScale(width[2], width[1], width[0])
            if self.glda:
                GL.glDrawArrays(GL.GL_QUADS, 24*i, 24)
            else:
                GL.glBegin(GL.GL_QUADS)
                for fi in range(6):
                    for v in range(4):
                        #print -LE, width
                        off = i*72+fi*12+v*3
                        qv = self.v[off:off+3]
                        tv = self.t[off:off+3]
                        GL.glTexCoord3f(*tv)
                        GL.glVertex3f(*qv)
                GL.glEnd()
            GL.glPopMatrix()
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        GL.glDisableClientState(GL.GL_TEXTURE_COORD_ARRAY)
        GL.glDisableClientState(GL.GL_TEXTURE_COORD_ARRAY_EXT)
        self.vertices.unbind()
        self.tvertices.unbind()
        GL.glDisable(GL.GL_CULL_FACE)

    def draw(self):

        self._draw_boxes()

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glUseProgram(self.fbo_program)

        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()

        GL.glActiveTexture(GL.GL_TEXTURE0)
        GL.glEnable(GL.GL_TEXTURE_2D)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.gl_state["fbo_tex"])

        GL.glActiveTexture(GL.GL_TEXTURE1)
        GL.glEnable(GL.GL_TEXTURE_1D)
        GL.glBindTexture(GL.GL_TEXTURE_1D, self.gl_state["cmap_tex"])

        GL.glUniform1i(self.fbo_uniform_locations["buffer"], 0)
        GL.glUniform1i(self.fbo_uniform_locations["colormap"], 1)

        GL.glColor3f(0.3, 0.5, 1.0)

        GL.glBegin(GL.GL_QUADS)
        GL.glTexCoord2i(0, 0)
        GL.glVertex3i(-1, -1, -1)

        GL.glTexCoord2i(0, 1)
        GL.glVertex3i( 1, -1, -1)

        GL.glTexCoord2i(1, 1)
        GL.glVertex3i( 1,  1, -1)

        GL.glTexCoord2i(1, 0)
        GL.glVertex3i(-1,  1, -1)
        GL.glEnd()

        GL.glUseProgram(0)

        GLUT.glutSwapBuffers()

    def move_slice(self, value):
        self.coord += value

    def rotate(self, axis, value):
        self.rotation[axis] += value/self.rfac

    def reset_view(self):   
        print "RESETTING"
        self.position = np.zeros(3, dtype='float')
        self.rotation = np.zeros(3, dtype='float') + 30
        self.position[2] = -2 # Offset backwards a bit

    def translate(self, axis, value):
        self.position[axis] += value/self.tfac

    def toggle_wireframe(self):
        self.wireframe = not self.wireframe
        print "Wireframe:", self.wireframe
        
    def toggle_glda(self):
        self.glda = not self.glda
        print "GLDA:", self.glda
        
    def _setup_keypress_handler(self):
        self.keypress_handler = ViewHandler3D(self)
        self.keypress_handler.dispatch_table.update(dict(
            y = (self.move_slice, ( 0.05,)),
            h = (self.move_slice, (-0.05,)),
            t = (self.toggle_wireframe, ()),
            u = (self.toggle_glda, ()),
            i = (self.recompile, ()),
            o = (self.reset_view, ())
            ))

# We override the standard PyOpenGL one because otherwise we can't validate
# multitexture shaders.
def compileProgram(*my_shaders, **kwargs):
    """Create a new program, attach my_shaders and validate
    
    my_shaders -- arbitrary number of my_shaders to attach to the 
        generated program.
    
    This convenience function is *not* standard OpenGL,
    but it does wind up being fairly useful for demos 
    and the like.  You may wish to copy it to your code 
    base to guard against PyOpenGL changes.
    
    Usage:
    
        shader = compileProgram( 
            compileShader( source, GL_VERTEX_SHADER ),
            compileShader( source2, GL_FRAGMENT_SHADER ),
        )
        glUseProgram( shader )
    
    Note:
        If (and only if) validation of the linked program 
        *passes* then the passed-in shader objects will be 
        deleted from the GL.
    
    returns GLuint shader program reference
    raises RuntimeError when a link/validation failure occurs
    """
    program = GL.glCreateProgram()
    for shader in my_shaders:
        GL.glAttachShader(program, shader)
    GL.glLinkProgram(program)
    # Validation has to occur *after* linking
    if 'callback' in kwargs: kwargs['callback'](program)
    GL.glValidateProgram( program )
    validation = GL.glGetProgramiv( program, GL.GL_VALIDATE_STATUS )
    if validation == GL.GL_FALSE:
        raise RuntimeError(
            """Validation failure (%s): %s"""%(
            validation,
            GL.glGetProgramInfoLog( program ),
        ))
    link_status = GL.glGetProgramiv( program, GL.GL_LINK_STATUS )
    if link_status == GL.GL_FALSE:
        raise RuntimeError(
            """Link failure (%s): %s"""%(
            link_status,
            GL.glGetProgramInfoLog( program ),
        ))
    for shader in my_shaders:
        GL.glDeleteShader(shader)
    return shaders.ShaderProgram( program )

if __name__ == "__main__":
    import yt.convenience, yt.frontends.enzo.api
    pf = yt.convenience.load(sys.argv[-1])
    print pf
    hv = HomogenizedVolume(pf = pf, fields=["Density"], log_fields=[True])
    hv.initialize_source()
    mip = MIPScene(hv)
    mip.run()
