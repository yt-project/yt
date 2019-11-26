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

import cyglfw3 as glfw
import numpy as np
import OpenGL.GL as GL
from .input_events import EventCollection, MouseRotation, JoystickAction
from .glfw_inputhook import InputHookGLFW
from ..image_writer import write_bitmap

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

class RenderingContext(object):
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
    should_quit = False
    image_widget = None
    draw = False
    def __init__(self, width=1024, height=1024, title="vol_render",
                 always_on_top = False, decorated = True, position = None,
                 visible = True):
        curdir = os.getcwd()
        glfw.Init()
        # glfw sometimes changes the current working directory, see
        # https://github.com/adamlwgriffiths/cyglfw3/issues/21
        os.chdir(curdir)
        self.offscreen = not visible
        glfw.WindowHint(glfw.VISIBLE, visible)
        glfw.WindowHint(glfw.CONTEXT_VERSION_MAJOR, 3)
        glfw.WindowHint(glfw.CONTEXT_VERSION_MINOR, 3)
        glfw.WindowHint(glfw.OPENGL_FORWARD_COMPAT, True)
        glfw.WindowHint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
        glfw.WindowHint(glfw.FLOATING, always_on_top)
        glfw.WindowHint(glfw.DECORATED, decorated)
        self.window = glfw.CreateWindow(width, height, title)
        if not self.window:
            glfw.Terminate()
            exit()
        if position is None:
            self.center_window()
        else:
            self.set_position(*position)

        glfw.MakeContextCurrent(self.window)
        GL.glClearColor(0.0, 0.0, 0.0, 0.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        glfw.SwapBuffers(self.window)
        glfw.PollEvents()

    def set_position(self, xpos, ypos):
        if xpos < 0 or ypos < 0:
            raise RuntimeError
        monitor = glfw.GetPrimaryMonitor()
        video_mode = glfw.GetVideoMode(monitor)
        win_width, win_height = glfw.GetWindowSize(self.window)
        if 0 < xpos < 1:
            # We're being fed relative coords.  We offset these for the window
            # center.
            xpos = max(xpos * video_mode.width - 0.5 * win_width, 0)
        if 0 < ypos < 1:
            # We're being fed relative coords.  We offset these for the window
            # center.
            ypos = max(ypos * video_mode.height - 0.5 * win_height, 0)
        print("Setting position", xpos, ypos)
        glfw.SetWindowPos(self.window, xpos, ypos)

    def center_window(self):
        self.set_position(0.5, 0.5)

    def setup_loop(self, scene, camera):
        scene.camera = camera
        camera.compute_matrices()
        print("Starting rendering...")
        callbacks = EventCollection(scene, camera)
        # register key callbacks defined in 
        # yt.visualization.volume_rendering.input_events
        callbacks.add_key_callback("close_window", "escape")
        callbacks.add_key_callback("zoomin", "w")
        callbacks.add_key_callback("zoomout", "s")
        callbacks.add_key_callback("closeup", "z")
        callbacks.add_key_callback("cmap_cycle", "c")
        callbacks.add_key_callback("reset", "r")
        callbacks.add_key_callback("camera_orto", "o")
        callbacks.add_key_callback("camera_proj", "p")
        callbacks.add_key_callback("shader_max", "1")
        callbacks.add_key_callback("shader_proj", "2")
        callbacks.add_key_callback("shader_test", "3")
        callbacks.add_key_callback("shader_lines", "4")
        callbacks.add_key_callback("print_limits", "g")
        callbacks.add_key_callback("print_help", "h")
        callbacks.add_key_callback("debug_buffer", "d")
        callbacks.add_key_callback("cmap_max_up", "right_bracket")
        callbacks.add_key_callback("cmap_max_down", "left_bracket")
        callbacks.add_key_callback("cmap_min_up", "semicolon")
        callbacks.add_key_callback("cmap_min_down", "apostrophe")
        callbacks.add_key_callback("cmap_autoscale", "a")
        callbacks.add_key_callback("cmap_toggle_log", "l")
        callbacks.add_key_callback("screenshot", "f")
        callbacks.add_key_callback("prev_component", "minus")
        callbacks.add_key_callback("next_component", "equal")
        mouse_callbacks = MouseRotation()
        callbacks.add_mouse_callback(mouse_callbacks.start_rotation,
            glfw.MOUSE_BUTTON_LEFT)
        callbacks.add_mouse_callback(mouse_callbacks.stop_rotation,
            glfw.MOUSE_BUTTON_LEFT, action="release")
        callbacks.add_framebuffer_callback("framebuffer_size")
        callbacks.add_render_callback(mouse_callbacks.do_rotation)
        joystick_callbacks = JoystickAction()
        #joystick_callbacks.calibrate()
        callbacks.add_render_callback(joystick_callbacks.check_axes)
        glfw.SetFramebufferSizeCallback(self.window,
            callbacks.framebuffer_call)
        glfw.SetKeyCallback(self.window, callbacks.key_call)
        glfw.SetMouseButtonCallback(self.window, callbacks.mouse_call)
        callbacks.draw = True
        return callbacks

    def start_shell(self, scene, camera, **kwargs): 
        from IPython.terminal.pt_inputhooks import register
        from IPython.terminal.embed import InteractiveShellEmbed
        callbacks = self.setup_loop(scene, camera)
        rl = self(scene, camera, callbacks)
        hook = InputHookGLFW(rl)
        register("glfw", hook)
        shell = InteractiveShellEmbed(**kwargs)
        shell.enable_gui('glfw')
        return shell

    def start_loop(self, scene, camera):
        callbacks = self.setup_loop(scene, camera)
        for i in self(scene, camera, callbacks):
            pass

    def __call__(self, scene, camera, callbacks):
        while not glfw.WindowShouldClose(self.window) or self.should_quit:
            callbacks(self.window)
            if callbacks.draw or self.draw:
                camera.compute_matrices()
                scene.camera = camera
                scene.render()
                glfw.SwapBuffers(self.window)
                callbacks.draw = self.draw = False
                if self.image_widget is not None:
                    self.image_widget.value = write_bitmap(
                            scene.image[:,:,:3], None)
            glfw.PollEvents()
            yield self
        glfw.Terminate()

