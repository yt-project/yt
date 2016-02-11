from collections import defaultdict, namedtuple
from yt.utilities.math_utils import \
    get_perspective_matrix, \
    get_orthographic_matrix
from yt.utilities.exceptions import YTInvalidShaderType
from functools import wraps
import OpenGL.GL as GL
import cyglfw3 as glfw
import numpy as np
import matplotlib.cm as cm
import random

event_registry = {}

GLFWEvent = namedtuple("GLFWEvent", ['window', 'key', 'scancode', 'action',
                       'mods', 'width', 'height'])

class EventCollection(object):
    def __init__(self, camera):
        self.key_callbacks = defaultdict(list)
        self.mouse_callbacks = defaultdict(list)
        self.framebuffer_callbacks = []
        self.render_events = []
        self.camera = camera
        self.draw = True

    def key_call(self, window, key, scancode, action, mods):
        draw = False
        event = GLFWEvent(window, key, scancode, action, mods, None, None)
        for f in self.key_callbacks[key, action, mods]:
            draw = f(self, event) or draw
        self.draw = self.draw or draw

    def mouse_call(self, window, key, action, mods):
        event = GLFWEvent(window, key, None, action, mods, None, None)
        draw = False
        for f in self.mouse_callbacks[key, action, mods]:
            draw = f(self, event) or draw
        self.draw = self.draw or draw

    def framebuffer_call(self, window, width, height):
        event = GLFWEvent(window, None, None, None, None, width, height)
        draw = False
        for f in self.framebuffer_callbacks:
            draw = f(self, event) or draw
        self.draw = self.draw or draw

    def __call__(self, window):
        event = GLFWEvent(window, None, None, None, None, None, None)
        draw = False
        for f in self.render_events:
            draw = f(self, event) or draw
        self.draw = self.draw or draw

    def add_render_callback(self, func):
        self.render_events.append(func)

    def add_key_callback(self, func, key, action = "press", mods = None):
        self._add_callback(self.key_callbacks, func, key, action, mods)

    def add_mouse_callback(self, func, key, action = "press", mods = None):
        self._add_callback(self.mouse_callbacks, func, key, action, mods)

    def add_framebuffer_callback(self, func):
        self.framebuffer_callbacks.append(func)

    def _add_callback(self, d, func, key, action, mods):
        if not callable(func):
            func = event_registry[func]
        if isinstance(key, str):
            key = getattr(glfw, "KEY_%s" % key.upper())
        if isinstance(action, str):
            action = getattr(glfw, action.upper())
        if not isinstance(mods, tuple):
            mods = (mods, )
        mod = 0
        for m in mods:
            if isinstance(m, str):
                m = getattr(glfw, "MOD_%s" % m.upper())
            elif m is None:
                m = 0
            mod |= m
        # We can allow for multiple
        d[key, action, mod].append(func)

def register_event(name):
    def _f(func):
        event_registry[name] = func
        return func
    return _f

@register_event("framebuffer_size")
def framebuffer_size_callback(event_coll, event):
    GL.glViewport(0, 0, event.width, event.height)
    event_coll.camera.aspect_ratio = float(event.width)/event.height
    return True

@register_event("close_window")
def close_window(event_coll, event):
    glfw.SetWindowShouldClose(event.window, True)

@register_event("zoomin")
def zoomin(event_coll, event):
    camera = event_coll.camera
    camera.position -= 0.05 * (camera.position - camera.focus) / \
                np.linalg.norm(camera.position - camera.focus)
    return True

@register_event("zoomout")
def zoomout(event_coll, event):
    camera = event_coll.camera
    camera.position += 0.05 * (camera.position - camera.focus) / \
        np.linalg.norm(camera.position - camera.focus)
    return True

@register_event("camera_orto")
def camera_orto(event_coll, event):
    camera = event_coll.camera
    if camera.proj_func == get_orthographic_matrix:
        return False
    camera.proj_func = get_orthographic_matrix
    camera.fov = np.tan(np.radians(camera.fov) / 2.0)
    return True

@register_event("camera_proj")
def camera_proj(event_coll, event):
    camera = event_coll.camera
    if camera.proj_func == get_perspective_matrix:
        return False
    camera.proj_func = get_perspective_matrix
    camera.fov = np.degrees(np.arctan(camera.fov) * 2.0)
    return True

@register_event("cmap_cycle")
def cmap_cycle(event_coll, event):
    cmap = ['algae', 'kamae', 'viridis', 'inferno', 'magma']
    cmap = cm.get_cmap(random.choice(cmap))
    event_coll.camera.cmap = np.array(cmap(np.linspace(0, 1, 256)),
        dtype=np.float32)
    event_coll.camera.cmap_new = True
    print("Setting colormap to {}".format(cmap.name))
    return True

@register_event("closeup")
def closeup(event_coll, event):
    event_coll.camera.position = (0.01, 0.01, 0.01)
    return True

@register_event("reset")
def reset(event_coll, event):
    event_coll.camera.position = (-1.0, -1.0, -1.0)
    return True

class MouseRotation(object):
    def __init__(self):
        self.start = None
        self.rotation = False

    def start_rotation(self, event_coll, event):
        start_screen = glfw.GetCursorPos(event.window) # Screen coordinates
        window_size = glfw.GetWindowSize(event.window)

        norm_x = -1.0 + 2.0 * start_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * start_screen[1] / window_size[1]
        self.start = (norm_x, norm_y)
        self.rotation = True
        return False

    def stop_rotation(self, event_coll, event):
        end_screen = glfw.GetCursorPos(event.window)
        window_size = glfw.GetWindowSize(event.window)

        norm_x = -1.0 + 2.0 * end_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * end_screen[1] / window_size[1]
        end = (norm_x, norm_y)

        event_coll.camera.update_orientation(
            self.start[0], self.start[1], end[0], end[1])
        self.rotation = False
        return True

    def do_rotation(self, event_coll, event):
        if not self.rotation: return False
        new_end_screen = glfw.GetCursorPos(event.window)
        window_size = glfw.GetWindowSize(event.window)

        norm_x = -1.0 + 2.0 * new_end_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * new_end_screen[1] / window_size[1]
        new_end = (norm_x, norm_y)

        event_coll.camera.update_orientation(
            self.start[0], self.start[1], new_end[0], new_end[1])
        self.start = new_end
        return True

