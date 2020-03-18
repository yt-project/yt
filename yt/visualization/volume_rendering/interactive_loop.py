# encoding: utf-8
"""
Event loop for Interactive Data Visualization

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

import pyglet
import numpy as np
import OpenGL.GL as GL
from .input_events import EventCollection, MouseRotation, JoystickAction
from .glfw_inputhook import InputHookGLFW
from ..image_writer import write_bitmap
try:
    from .simple_gui import SimpleGUI
except ImportError:
    pass

from yt import write_bitmap

class EGLRenderingContext(object):
    '''Rendering context using EGL (experimental)

    Parameters
    ----------
    width : int, optional
        The width of the off-screen buffer window.  For performance reasons it
        is recommended to use values that are natural powers of 2.

    height : int, optional
        The height of the off-screen buffer window.  For performance reasons it
        it is recommended to use values that are natural powers of 2.

    '''

    def __init__(self, width=1024, height=1024):
        from OpenGL import EGL
        self.EGL = EGL
        self.display = EGL.eglGetDisplay(EGL.EGL_DEFAULT_DISPLAY)
        major = np.zeros(1, "i4")
        minor = np.zeros(1, "i4")
        EGL.eglInitialize(self.display, major, minor)
        num_configs = np.zeros(1, "i4")
        config = EGL.EGLConfig()
        # Now we create our necessary bits.
        config_attribs = np.array([
          EGL.EGL_SURFACE_TYPE, EGL.EGL_PBUFFER_BIT,
          EGL.EGL_BLUE_SIZE, 8,
          EGL.EGL_GREEN_SIZE, 8,
          EGL.EGL_RED_SIZE, 8,
          EGL.EGL_DEPTH_SIZE, 8,
          EGL.EGL_RENDERABLE_TYPE,
          EGL.EGL_OPENGL_BIT,
          EGL.EGL_NONE,
        ], dtype="i4")
        self.config = EGL.eglChooseConfig(self.display, config_attribs, config, 1,
            num_configs)

        pbuffer_attribs = np.array([
          EGL.EGL_WIDTH, width,
          EGL.EGL_HEIGHT, height,
          EGL.EGL_NONE
        ], dtype="i4")
        self.surface = EGL.eglCreatePbufferSurface(self.display, self.config,
            pbuffer_attribs)
        EGL.eglBindAPI(EGL.EGL_OPENGL_API)
        
        self.context = EGL.eglCreateContext(self.display, self.config,
            EGL.EGL_NO_CONTEXT, None)

        EGL.eglMakeCurrent(self.display, self.surface, self.surface,
            self.context)

        GL.glClearColor(0.0, 0.0, 0.0, 0.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

    def setup_loop(self, scene, camera):
        scene.set_camera(camera)
        scene.update_minmax()
        camera.compute_matrices()
        callbacks = EventCollection(scene, camera)
        callbacks.draw = True
        return callbacks

    def start_loop(self, scene, camera):
        callbacks = self.setup_loop(scene, camera)
        for i in self(scene, camera, callbacks):
            pass

    def __call__(self, scene, camera, callbacks):
        camera.compute_matrices()
        scene.set_camera(camera)
        scene.render()
        arr = scene._retrieve_framebuffer()
        write_bitmap(arr, "test.png")

class RenderingContext(pyglet.window.Window):
    '''Basic rendering context for IDV using GLFW3, that handles the main window even loop
    
    Parameters
    ----------
    width : int, optional
        The width of the Interactive Data Visualization window.  For
        performance reasons it is recommended to use values that are natural
        powers of 2.
    height : int, optional
        The height of the Interactive Data Visualization window.  For
        performance reasons it is recommended to use values that are natural
        powers of 2.
    title : str, optional
        The title of the Interactive Data Visualization window. 
    always_on_top : bool, optional
        Should this window be created such that it is always on top of other
        windows? (Default: False)
    decorated : bool, optional
        Does the window have operating system widgets (minimize, maximize
        close), or is it a bare context? (Default: True)
    position : tuple of ints, optional
        What position should the window be moved to? (Upper left)  If not
        specified, default to center.
    '''
    image_widget = None
    def __init__(self, width=1024, height=1024, title="vol_render",
                 always_on_top = False, decorated = True, position = None,
                 visible = True, scene = None):
        self.offscreen = not visible
        config = pyglet.gl.Config(major_version = 3, minor_version = 3,
                                  forward_compat = True, double_buffer = True)
        super(RenderingContext, self).__init__(width, height,
                                           config = config,
                                           visible = visible,
                                           caption = title)
        if position is None:
            self.center_window()
        else:
            #self.set_position(*position)
            self.set_location(*position)

        self.scene = scene
        self.label = pyglet.text.Label('Hello, yt',
                                  font_name='Times New Roman',
                                  font_size=36,
                                  x=self.width//2, y=self.height//2,
                                  anchor_x='center', anchor_y='center')

    def on_draw(self):
        self.clear()
        if self.scene is not None:
            self.scene.render()
            if self.image_widget is not None:
                self.image_widget.value = write_bitmap(
                        self.scene.image[:,:,:3], None)
        self.label.draw()

    def set_position(self, xpos, ypos):
        if xpos < 0 or ypos < 0:
            raise RuntimeError
        max_width = self.screen.width
        max_height = self.screen.height
        win_width, win_height = self.width, self.height
        if 0 < xpos < 1:
            # We're being fed relative coords.  We offset these for the window
            # center.
            xpos = max(xpos * max_width - 0.5 * win_width, 0)
        if 0 < ypos < 1:
            # We're being fed relative coords.  We offset these for the window
            # center.
            ypos = max(ypos * max_height - 0.5 * win_height, 0)
        print("Setting position", xpos, ypos)
        self.set_location(int(xpos), int(ypos))

    def center_window(self):
        self.set_position(0.5, 0.5)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        print(x, y, dx, dy, buttons, modifiers)
