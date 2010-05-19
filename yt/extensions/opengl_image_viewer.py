"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation:  UCSD
License:
  Copyright (C) 2010 Matthew Turk  All Rights Reserved.

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
import Image
import glob
import numpy as na
import time

ESCAPE = '\033'

class MultiImageDisplayScene(object):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH)

    def __init__(self):
        self._frames = []
        self._current = -1
        self.init_glut(512, 512)
        self.init_opengl(512, 512)

    def add_image(self, obj):
        self._frames.append(obj)

    def init_glut(self, width, height):
        GLUT.glutInit([]) # drop sys.argv
        GLUT.glutInitDisplayMode(self._display_mode)
        GLUT.glutInitWindowSize(width, height)
        GLUT.glutInitWindowPosition(0, 0)
        window = GLUT.glutCreateWindow("Image Display")
        GLUT.glutDisplayFunc(self.draw)
        GLUT.glutIdleFunc(self.draw)
        GLUT.glutKeyboardFunc(self.keypress_handler)

    def run(self):
        GLUT.glutMainLoop()

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

class StereoMultiImageDisplayScene(MultiImageDisplayScene):
    _display_mode = (GLUT.GLUT_RGBA | GLUT.GLUT_DOUBLE | GLUT.GLUT_DEPTH |
                     GLUT.GLUT_STEREO)

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
        buffer = na.array(Image.open(fn))
        print "Uploading buffer", buffer.min(), buffer.max(), buffer.shape, buffer.dtype
        obj = cls(tex_unit)
        obj.upload_image(buffer)
        return obj

class StereoImagePair(FlatImage):
    def __init__(self, tex_unit = GL.GL_TEXTURE0):
        self.tex_unit = tex_unit
        self.left_image = FlatImage(tex_unit)
        self.right_image = FlatImage(tex_unit)

    def draw(self):
        # Look left
        GL.glDrawBuffer(GL.GL_BACK_LEFT)
        self.left_image.draw()
        # Look right
        GL.glDrawBuffer(GL.GL_BACK_RIGHT)
        self.right_image.draw()

    def upload_images(self, buffer_left, buffer_right):
        self.left_image.upload_image(buffer_left)
        self.right_image.upload_image(buffer_right)

    @classmethod
    def from_image_files(cls, left_fn, right_fn, tex_unit = GL.GL_TEXTURE0):
        print "Uploading pairs from %s and %s" % (left_fn, right_fn)
        left_buffer = na.array(Image.open(left_fn))
        right_buffer = na.array(Image.open(right_fn))
        obj = cls(tex_unit)
        obj.left_image.upload_image(left_buffer)
        obj.right_image.upload_image(right_buffer)
        return obj

if __name__ == "__main__":
    fn_list = glob.glob("frames/*.png")

    main_scene = MultiImageDisplayScene()
    for fn in sorted(fn_list):
        main_scene.add_image(FlatImage.from_image_file(fn))
    main_scene._current = 0
    main_scene.run()
