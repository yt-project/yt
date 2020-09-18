# encoding: utf-8
"""
Event loop for Interactive Data Visualization

"""

import numpy as np
import pyglet
from OpenGL import GL

from ..image_writer import write_bitmap
from .input_events import EventCollection

try:
    from .simple_gui import SimpleGUI
except ImportError:
    pass


class EGLRenderingContext:
    """Rendering context using EGL (experimental)

    Parameters
    ----------
    width : int, optional
        The width of the off-screen buffer window.  For performance reasons it
        is recommended to use values that are natural powers of 2.

    height : int, optional
        The height of the off-screen buffer window.  For performance reasons it
        it is recommended to use values that are natural powers of 2.

    """

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
        config_attribs = np.array(
            [
                EGL.EGL_SURFACE_TYPE,
                EGL.EGL_PBUFFER_BIT,
                EGL.EGL_BLUE_SIZE,
                8,
                EGL.EGL_GREEN_SIZE,
                8,
                EGL.EGL_RED_SIZE,
                8,
                EGL.EGL_DEPTH_SIZE,
                8,
                EGL.EGL_RENDERABLE_TYPE,
                EGL.EGL_OPENGL_BIT,
                EGL.EGL_NONE,
            ],
            dtype="i4",
        )
        self.config = EGL.eglChooseConfig(
            self.display, config_attribs, config, 1, num_configs
        )

        pbuffer_attribs = np.array(
            [EGL.EGL_WIDTH, width, EGL.EGL_HEIGHT, height, EGL.EGL_NONE], dtype="i4"
        )
        self.surface = EGL.eglCreatePbufferSurface(
            self.display, self.config, pbuffer_attribs
        )
        EGL.eglBindAPI(EGL.EGL_OPENGL_API)

        self.context = EGL.eglCreateContext(
            self.display, self.config, EGL.EGL_NO_CONTEXT, None
        )

        EGL.eglMakeCurrent(self.display, self.surface, self.surface, self.context)

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
        self.setup_loop(scene, camera)

    def __call__(self, scene, camera, callbacks):
        camera.compute_matrices()
        scene.set_camera(camera)
        scene.render()
        arr = scene._retrieve_framebuffer()
        write_bitmap(arr, "test.png")


class RenderingContext(pyglet.window.Window):
    """Basic rendering context for IDV using GLFW3, that handles the main window even loop

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
    """

    image_widget = None
    _do_update = True

    def __init__(
        self,
        width=1024,
        height=1024,
        title="vol_render",
        always_on_top=False,
        decorated=True,
        position=None,
        visible=True,
        gui=False,
    ):
        self.offscreen = not visible
        config = pyglet.gl.Config(
            major_version=3, minor_version=3, forward_compat=True, double_buffer=True
        )
        super(RenderingContext, self).__init__(
            width, height, config=config, visible=visible, caption=title
        )
        if position is None:
            self.center_window()
        else:
            # self.set_position(*position)
            self.set_location(*position)

        if gui:
            gui = SimpleGUI(self)
        self.gui = gui

    def on_draw(self):
        self.switch_to()
        if self._do_update and self.scene:
            self._do_update = False
            self.clear()
            if self.scene is not None:
                # This should check if the scene is actually dirty, and only
                # redraw the FB if it's not; this might need another flag
                self.scene.render()
                if self.image_widget is not None:
                    self.image_widget.value = write_bitmap(
                        self.scene.image[:, :, :3], None
                    )
        if self.gui:
            self.switch_to()
            self.gui.render(self.scene)

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

    def on_mouse_press(self, x, y, button, modifiers):
        self._do_update = True

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if self.gui and self.gui.mouse_event_handled:
            self._do_update = True
            return
        start_x = -1.0 + 2.0 * x / self.width
        end_x = -1.0 + 2.0 * (x - dx) / self.width
        start_y = -1.0 + 2.0 * y / self.height
        end_y = -1.0 + 2.0 * (y + dy) / self.height

        self.scene.camera.update_orientation(start_x, start_y, end_x, end_y)
        self._do_update = True

    def on_mouse_scroll(self, x, y, scroll_x, scroll_y):
        # captures mouse scrolling as zoom in/out
        if self.gui and self.gui.mouse_event_handled:
            self._do_update = True
            return

        camera = self.scene.camera  # current camera
        dPos = (
            0.1
            * (camera.position - camera.focus)
            / np.linalg.norm(camera.position - camera.focus)
        )

        # wheel scroll comes in the scroll_y parameter with a value +/- 1
        # +1 when scrolling "down", -1 when scrolling "up", so
        # flip it so scrolling "down" zooms out:
        zoom_inout = -1 * scroll_y
        self.scene.camera.offsetPosition(zoom_inout * dPos)
        self._do_update = True

    def on_key_press(self, symbol, modifiers):
        # skeleton for capturing key presses!

        # potential navigation keys
        if symbol in [pyglet.window.key.LEFT, pyglet.window.key.A]:
            pass
        if symbol in [pyglet.window.key.RIGHT, pyglet.window.key.D]:
            pass
        if symbol in [pyglet.window.key.UP, pyglet.window.key.W]:
            pass
        if symbol in [pyglet.window.key.DOWN, pyglet.window.key.S]:
            pass
