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

import sys
import OpenGL.GL as GL
import OpenGL.GLUT as GLUT
import OpenGL.GLU as GLU
import OpenGL.GL.shaders as shaders
from OpenGL.arrays import vbo, ArrayDatatype
import Image
import glob
import numpy as np
import time

ESCAPE = '\033'

class ViewHandler3D(object):
    def __init__(self, scene):
        # We 
        self.scene = scene
        self.dispatch_table = dict(
            q = (scene.translate, (1,  1.0)),
            e = (scene.translate, (1, -1.0)),
            w = (scene.translate, (2,  1.0)),
            s = (scene.translate, (2, -1.0)),
            a = (scene.translate, (0,  1.0)),
            d = (scene.translate, (0, -1.0)),

            Q = (scene.rotate, (1,  1.0)),
            E = (scene.rotate, (1, -1.0)),
            W = (scene.rotate, (2,  1.0)),
            S = (scene.rotate, (2, -1.0)),
            A = (scene.rotate, (0,  1.0)),
            D = (scene.rotate, (0, -1.0)),

            ESCAPE = (sys.exit, (0,))
        )

    def __call__(self, *args):
        # We set up our standard handlers, and then anything additional can get
        # called if none of our dispatch mechanisms work.
        if args[0] in self.dispatch_table:
            func, args = self.dispatch_table[args[0]]
            func(*args)
        # always draw when handling a keypress, even if it's one time too many
        self.scene.draw() 

class GenericGLUTScene(object):
    
    def __init__(self, width, height):
        self.init_glut(width, height)
        self.init_opengl(width, height)

    def init_glut(self, width, height):
        GLUT.glutInit([]) # drop sys.argv
        GLUT.glutInitDisplayMode(self._display_mode)
        GLUT.glutInitWindowSize(width, height)
        GLUT.glutInitWindowPosition(0, 0)
        self.window = GLUT.glutCreateWindow(self._title)
        GLUT.glutDisplayFunc(self.draw)
        #GLUT.glutIdleFunc(self.draw)
        GLUT.glutKeyboardFunc(self.keypress_handler)

    def run(self):
        GLUT.glutMainLoop()

class MultiImageDisplayScene(object):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH)
    _title = "Image Display"
    def __init__(self):
        GenericGLUTScene(512, 512)
        self._frames = []
        self._current = -1

    def add_image(self, obj):
        self._frames.append(obj)

    def init_opengl(self, width, height):
        GL.glClearColor(0.0, 0.0, 0.0, 0.0)
        GL.glClearDepth(1.0)
        GL.glDepthFunc(GL.GL_LESS)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glShadeModel(GL.GL_SMOOTH)
        
    def draw(self):

        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        if self._current >= 0:
            self._draw_current()

        GLUT.glutSwapBuffers()

    def _draw_current(self):
        self._frames[self._current].draw()

    def keypress_handler(self, *args):
        if args[0] == ESCAPE:
            sys.exit()
        elif len(self._frames) == 0:
            # The rest only operate if we have multiple frames
            return
        elif args[0] == 'n':
            self._current = min(self._current + 1, len(self._frames))
        elif args[0] == 'N':
            while self._current < len(self._frames) - 1:
                self._current += 1
                self.draw()
        elif args[0] == 'p':
            self._current = max(self._current - 1, 0)
        elif args[0] == 'P':
            while self._current > 0:
                self._current -= 1
                self.draw()
        elif args[0] == 't':
            for i in xrange(15):
                time.sleep(0.05)
                self.draw()
                self._current += 1
                time.sleep(0.05)
                self.draw()
                self._current -= 1
        self.draw() # Once more for good measure

class StereoMultiImageDisplayScene(MultiImageDisplayScene):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH |
                     GLUT.GLUT_STEREO)

    def draw(self):
        GL.glDrawBuffer(GL.GL_BACK_LEFT)

        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        if self._current >= 0:
            self._draw_current_left()

        GLUT.glutSwapBuffers()

        GL.glDrawBuffer(GL.GL_BACK_RIGHT)

        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        if self._current >= 0:
            self._draw_current_right()

        GLUT.glutSwapBuffers()

    def _draw_current_left(self):
        self._frames[self._current].draw_left()

    def _draw_current_right(self):
        self._frames[self._current].draw_right()

class FlatImage(object):
    def __init__(self, tex_unit = GL.GL_TEXTURE0):
        self.tex_unit = tex_unit
        GL.glActiveTexture(self.tex_unit)
        self._id = GL.glGenTextures(1)

    def draw(self):
        GL.glActiveTexture(self.tex_unit)
        GL.glEnable(GL.GL_TEXTURE_2D)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self._id)
        GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_DECAL)

        GL.glColor3f(0.3, 0.5, 1.0)

        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3i(-1, -1, -1)
        GL.glTexCoord2i(0, 0)

        GL.glVertex3i( 1, -1, -1)
        GL.glTexCoord2i(0, 1)

        GL.glVertex3i( 1,  1, -1)
        GL.glTexCoord2i(1, 1)

        GL.glVertex3i(-1,  1, -1)
        GL.glTexCoord2i(1, 0)
        GL.glEnd()

    def upload_image(self, buffer):
        ix, iy, nc = buffer.shape
        if nc != 4: raise RuntimeError
        GL.glActiveTexture(self.tex_unit)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self._id)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, ix, iy, 0, GL.GL_RGBA,
                        GL.GL_UNSIGNED_BYTE, buffer.tostring())
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)

    @classmethod
    def from_image_file(cls, fn, tex_unit = GL.GL_TEXTURE0):
        buffer = np.array(Image.open(fn))
        print "Uploading buffer", buffer.min(), buffer.max(), buffer.shape, buffer.dtype
        obj = cls(tex_unit)
        obj.upload_image(buffer)
        return obj

class StereoImagePair(FlatImage):
    def __init__(self, tex_unit = GL.GL_TEXTURE0):
        self.tex_unit = tex_unit
        self.left_image = FlatImage(tex_unit)
        self.right_image = FlatImage(tex_unit)

    def draw_left(self):
        self.left_image.draw()

    def draw_right(self):
        self.right_image.draw()

    def upload_images(self, buffer_left, buffer_right):
        self.left_image.upload_image(buffer_left)
        self.right_image.upload_image(buffer_right)

    @classmethod
    def from_image_files(cls, left_fn, right_fn, tex_unit = GL.GL_TEXTURE0):
        print "Uploading pairs from %s and %s" % (left_fn, right_fn)
        left_buffer = np.array(Image.open(left_fn))
        right_buffer = np.array(Image.open(right_fn))
        obj = cls(tex_unit)
        obj.left_image.upload_image(left_buffer)
        obj.right_image.upload_image(right_buffer)
        return obj

_verts = ( (0,0,0), (1,0,0), (1,1,0), (0,1,0),
           (0,0,0), (1,0,0), (1,0,1), (0,0,1),
           (0,0,0), (0,0,1), (0,1,1), (0,1,0),
           (0,1,0), (1,1,0), (1,1,1), (0,1,1),
           (1,1,0), (1,0,0), (1,0,1), (1,1,1),
           (0,0,1), (0,1,1), (1,1,1), (1,0,1) )

class GridObject3DScene(GenericGLUTScene):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH)
    _title = "Grids"

    def _get_grid_vertices(self, offset):
        DLE, DRE = self.pf.domain_left_edge, pf.domain_right_edge
        DW = DRE - DLE
        k = 0
        for g in self.pf.h.grids:
            vs = ((g.LeftEdge-DLE)/DW, (g.RightEdge-DLE)/DW)
            for vert in _verts:
                for i,v in enumerate(vert):
                    yield vs[v][i] - offset
                    k += 1

    def __init__(self, pf, offset = 0.5):
        self.pf = pf
        GenericGLUTScene.__init__(self, 800, 800)

        num = len(pf.h.grids) * 6 * 4
        self.v = np.fromiter(self._get_grid_vertices(offset),
                             dtype = 'float32', count = num * 3)

        self.vertices = vbo.VBO(self.v)
        self.ng = len(pf.h.grids)
        self.ox = self.oy = self.rx = self.ry = self.rz = 0
        self.oz = -4

    def init_opengl(self, width, height):
        # One-time GL setup
        GL.glClearColor(1, 1, 1, 1)
        GL.glColor3f(1, 0, 0)
        GL.glEnable(GL.GL_DEPTH_TEST)
        #glEnable(GL_CULL_FACE)

        # Uncomment this line for a wireframe view
        GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE)

        # Simple light setup.  On Windows GL_LIGHT0 is enabled by default,
        # but this is not the case on Linux or Mac, so remember to always 
        # include it.
        GL.glEnable(GL.GL_LIGHTING)

        def vec(*args):
            return (GL.GLfloat * len(args))(*args)

        GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, vec(.5, .5, 1, 0))
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_SPECULAR, vec(.5, .5, 1, 1))
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, vec(1, 1, 1, 1))
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, vec(1, 0, .5, 0))
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE, vec(.5, .5, .5, 1))
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, vec(1, 1, 1, 1))

        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT_AND_DIFFUSE, vec(0.5, 0, 0.3, 1))
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, vec(1, 1, 1, 1))
        GL.glMaterialf(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, 50)

        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GLU.gluPerspective(60., width / float(height), 1e-3, 10.)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        
    def draw(self):

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glLoadIdentity()
        GL.glTranslatef(self.ox, self.oy, self.oz)
        GL.glRotatef(self.rx, 0, 0, 1)
        GL.glRotatef(self.ry, 0, 1, 0)
        GL.glRotatef(self.rz, 1, 0, 0)

        self.vertices.bind()
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glVertexPointer( 3, GL.GL_FLOAT, 0, self.vertices)
        GL.glDrawArrays(GL.GL_QUADS, 0, 4*6*self.ng)
        
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        self.vertices.unbind()
        GLUT.glutSwapBuffers()
        
    def keypress_handler(self, *args):
        tfac = 25.0
        rfac = 0.5
        if args[0] == ESCAPE:
            sys.exit()
        elif args[0] == 'a':
            self.ox += 1.0/tfac
        elif args[0] == 'd':
            self.ox -= 1.0/tfac
        elif args[0] == 's':
            self.oz -= 1.0/tfac
        elif args[0] == 'w':
            self.oz += 1.0/tfac
        elif args[0] == 'q':
            self.oy -= 1.0/tfac
        elif args[0] == 'e':
            self.oy += 1.0/tfac
        # Now, rotations
        elif args[0] == 'A':
            self.rx -= 1.0/rfac
        elif args[0] == 'D':
            self.rx += 1.0/rfac
        elif args[0] == 'S':
            self.rz -= 1.0/rfac
        elif args[0] == 'W':
            self.rz += 1.0/rfac
        elif args[0] == 'Q':
            self.ry -= 1.0/rfac
        elif args[0] == 'E':
            self.ry += 1.0/rfac
        self.draw()

class GridSlice3DScene(GenericGLUTScene):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH)
    _title = "Grids"

    def _get_grid_vertices(self, offset):
        for g in self.pf.h.grids:
            vs = (g.LeftEdge, g.RightEdge)
            for vert in _verts:
                for i,v in enumerate(vert):
                    yield vs[v][i] - offset

    def _setup_grids(self):
        self._grid_textures = {}
        for g in self.pf.h.grids:
            self._upload_grid_textures(g)

    def _upload_grid_textures(self, grid):
        ix, iy, iz = grid.ActiveDimensions

        GL.glActiveTexture(GL.GL_TEXTURE0)
        id_field = GL.glGenTextures(1)
        upload = np.log10(grid["Density"].astype("float32")).copy()
        self.mi = min(upload.min(), self.mi)
        self.ma = max(upload.max(), self.ma)
        #upload = (255*(upload - -31.0) / (-25.0 - -31.0)).astype("uint8")
        
        GL.glBindTexture(GL.GL_TEXTURE_3D, id_field)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
        GL.glTexImage3D(GL.GL_TEXTURE_3D, 0, GL.GL_LUMINANCE32F_ARB, iz, iy, ix, 0,
                        GL.GL_LUMINANCE, GL.GL_FLOAT, upload)

        GL.glActiveTexture(GL.GL_TEXTURE1)
        id_mask  = GL.glGenTextures(1)
        upload = grid.child_mask.astype("float32").copy()

        GL.glBindTexture(GL.GL_TEXTURE_3D, id_mask)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
        GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
        GL.glTexImage3D(GL.GL_TEXTURE_3D, 0, GL.GL_LUMINANCE, iz, iy, ix, 0,
                        GL.GL_LUMINANCE, GL.GL_FLOAT, upload)

        self._grid_textures[grid.id] = (id_field, id_mask)

        print "Uploaded", grid.id

    def __init__(self, pf, offset = 0.5):
        self.offset = offset
        self.mi, self.ma = 1e30, -1e30
        self.pf = pf
        self.coord = 0.0
        self.tfac = 10.0
        self.rfac = 0.5
        self._setup_keypress_handler()
        GenericGLUTScene.__init__(self, 800, 800)

        num = len(pf.h.grids) * 6 * 4
        self.v = np.fromiter(self._get_grid_vertices(offset),
                             dtype = 'float32', count = num * 3)

        self.vertices = vbo.VBO(self.v)
        self.ng = len(pf.h.grids)
        self.position = np.zeros(3, dtype='float')
        self.rotation = np.zeros(3, dtype='float')
        self.position[2] = -2 # Offset backwards a bit

        self._setup_grids()

    def init_opengl(self, width, height):
        # One-time GL setup
        GL.glClearColor(1, 1, 1, 1)
        GL.glColor3f(1, 0, 0)
        GL.glEnable(GL.GL_DEPTH_TEST)

        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GLU.gluPerspective(60., width / float(height), 1e-3, 10.)
        GL.glMatrixMode(GL.GL_MODELVIEW)

        # Now we compile our shaders

        self.program = shaders.compileProgram(
            shaders.compileShader('''
                void main() {
                    gl_TexCoord[0]=gl_TextureMatrix[0] * gl_MultiTexCoord0;
                    gl_TexCoord[1]=gl_TextureMatrix[1] * gl_MultiTexCoord1;
                    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
                }
            ''',GL.GL_VERTEX_SHADER),
            shaders.compileShader('''
                uniform float ma;
                uniform float mi;
                uniform sampler3D field;
                uniform sampler3D mask;

                void main() {
                    vec3 pos;
                    float val;

                    pos = vec3(gl_TexCoord[1].xyz);
                    val = texture3D( mask, pos )[0];
                    if(val == 0.0) discard;

                    pos = vec3(gl_TexCoord[0].xyz);
                    val = texture3D( field, pos )[0];

                    float color = (val - mi) / (ma - mi);
                
                    //gl_FragColor = vec4(pos.x, pos.y, pos.z, 1.0);
                    gl_FragColor = vec4(color, color, color, 1.0);
                }
        ''',GL.GL_FRAGMENT_SHADER),)

        self.uniform_locations = dict( (
                (v, GL.glGetUniformLocation(self.program, v)) for v in
                ['mi','ma','field','mask']
            ) )
            
    def draw(self):

        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE)

        GL.glLoadIdentity()
        GL.glTranslatef(*self.position)
        GL.glRotatef(self.rotation[0], 0, 0, 1)
        GL.glRotatef(self.rotation[1], 0, 1, 0)
        GL.glRotatef(self.rotation[2], 1, 0, 0)

        self.vertices.bind()
        GL.glColor3f(0.0, 0.0, 0.0)
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glVertexPointer( 3, GL.GL_FLOAT, 0, self.vertices)
        GL.glDrawArrays(GL.GL_QUADS, 0, 4*6*self.ng)
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        self.vertices.unbind()

        # Now, we just want to draw 
        GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)

        # We can just draw a single quad for now
        GL.glUseProgram(self.program)
        GL.glUniform1i(self.uniform_locations["field"], 0)
        GL.glUniform1i(self.uniform_locations["mask"], 1)
        GL.glUniform1f(self.uniform_locations["mi"], self.mi)
        GL.glUniform1f(self.uniform_locations["ma"], self.ma)
        GL.glEnable(GL.GL_TEXTURE_3D)

        t0, t1 = 0.0, 1.0
        for g in self.pf.h.find_slice_grids(self.coord + 0.5, 1)[0]:
            LE = g.LeftEdge - self.offset
            RE = g.RightEdge - self.offset
            off = (self.coord - LE[1]) / (RE[1] - LE[1])

            GL.glActiveTexture(GL.GL_TEXTURE0)
            GL.glBindTexture(GL.GL_TEXTURE_3D, self._grid_textures[g.id][0])
            GL.glEnable(GL.GL_TEXTURE_3D)

            GL.glActiveTexture(GL.GL_TEXTURE1)
            GL.glBindTexture(GL.GL_TEXTURE_3D, self._grid_textures[g.id][1])
            GL.glEnable(GL.GL_TEXTURE_3D)

            GL.glBegin(GL.GL_QUADS)

            GL.glMultiTexCoord3f(GL.GL_TEXTURE0, t0, off, t0)
            GL.glMultiTexCoord3f(GL.GL_TEXTURE1, t0, off, t0)
            GL.glVertex3f(LE[0], self.coord, LE[2])

            GL.glMultiTexCoord3f(GL.GL_TEXTURE0, t0, off, t1)
            GL.glMultiTexCoord3f(GL.GL_TEXTURE1, t0, off, t1)
            GL.glVertex3f(RE[0], self.coord, LE[2])

            GL.glMultiTexCoord3f(GL.GL_TEXTURE0, t1, off, t1)
            GL.glMultiTexCoord3f(GL.GL_TEXTURE1, t1, off, t1)
            GL.glVertex3f(RE[0], self.coord, RE[2])

            GL.glMultiTexCoord3f(GL.GL_TEXTURE0, t1, off, t0)
            GL.glMultiTexCoord3f(GL.GL_TEXTURE1, t1, off, t0)
            GL.glVertex3f(LE[0], self.coord, RE[2])

            GL.glEnd()

        GL.glUseProgram(0)

        GLUT.glutSwapBuffers()

    def move_slice(self, value):
        self.coord += value

    def rotate(self, axis, value):
        self.rotation[axis] += value/self.rfac

    def translate(self, axis, value):
        self.position[axis] += value/self.tfac
        
    def _setup_keypress_handler(self):
        self.keypress_handler = ViewHandler3D(self)
        self.keypress_handler.dispatch_table.update(dict(
            y = (self.move_slice, ( 0.05,)),
            h = (self.move_slice, (-0.05,))
            ))

if __name__ == "__main__":
    if sys.argv[-2] == '-g':
        import yt.mods
        pf = yt.mods.load(sys.argv[-1])
        main_scene = GridObject3DScene(pf)
    elif sys.argv[-2] == '-s':
        import yt.mods
        pf = yt.mods.load(sys.argv[-1])
        main_scene = GridSlice3DScene(pf)
    else:
        fn_list = glob.glob("frames/*.png")

        main_scene = MultiImageDisplayScene()
        for fn in sorted(fn_list):
            main_scene.add_image(FlatImage.from_image_file(fn))
        main_scene._current = 0
    main_scene.run()
