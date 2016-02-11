from collections import defaultdict
from functools import wraps
import OpenGL.GL as GL
import cyglfw3 as glfw
import numpy as np
import matplotlib.cm as cm
import random

event_registry = {}

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
        for f in self.key_callbacks[key, action, mods]:
            draw = f(self.camera, window, key, scancode, action, mods) or draw
        self.draw = self.draw or draw

    def mouse_call(self, window, key, action, mods):
        draw = False
        for f in self.mouse_callbacks[key, action, mods]:
            draw = f(self.camera, window, key, action, mods) or draw
        self.draw = self.draw or draw

    def framebuffer_call(self, window, width, height):
        draw = False
        for f in self.framebuffer_callbacks:
            draw = f(self.camera, window, width, height) or draw
        self.draw = self.draw or draw

    def __call__(self, window):
        draw = False
        for f in self.render_events:
            draw = f(self.camera, window) or draw
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
def framebuffer_size_callback(camera, window, width, height):
    GL.glViewport(0, 0, width, height)
    camera.aspect_ratio = float(width)/height
    return True

@register_event("close_window")
def close_window(camera, window, key, scancode, action, mods):
    glfw.SetWindowShouldClose(window, True)

@register_event("zoomin")
def zoomin(camera, window, key, scancode, action, mods):
    camera.position -= 0.05 * (camera.position - camera.focus) / \
                np.linalg.norm(camera.position - camera.focus)
    print camera.position, camera.focus
    return True

@register_event("zoomout")
def zoomout(camera, window, key, scancode, action, mods):
    camera.position += 0.05 * (camera.position - camera.focus) / \
        np.linalg.norm(camera.position - camera.focus)
    return True

@register_event("cmap_cycle")
def cmap_cycle(camera, window, key, scancode, action, mods):
    cmap = ['algae', 'kamae', 'viridis', 'inferno', 'magma']
    cmap = cm.get_cmap(random.choice(cmap))
    camera.cmap = np.array(cmap(np.linspace(0, 1, 256)), dtype=np.float32)
    camera.cmap_new = True
    print("Setting colormap to {}".format(cmap.name))
    return True

@register_event("closeup")
def closeup(camera, window, key, scancode, action, mods):
    camera.position = (0.01, 0.01, 0.01)
    return True

@register_event("reset")
def reset(camera, window, key, scancode, action, mods):
    camera.position = (-1.0, -1.0, -1.0)
    return True

@register_event("printit")
def printit(*args):
    print args

class MouseRotation(object):
    def __init__(self):
        self.start = None
        self.rotation = False

    def start_rotation(self, camera, window, key, action, mods):
        start_screen = glfw.GetCursorPos(window) # Screen coordinates
        window_size = glfw.GetWindowSize(window)

        norm_x = -1.0 + 2.0 * start_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * start_screen[1] / window_size[1]
        self.start = (norm_x, norm_y)
        self.rotation = True
        return False

    def stop_rotation(self, camera, window, key, action, mods):
        end_screen = glfw.GetCursorPos(window)
        window_size = glfw.GetWindowSize(window)

        norm_x = -1.0 + 2.0 * end_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * end_screen[1] / window_size[1]
        end = (norm_x, norm_y)

        camera.update_orientation(self.start[0],
                                  self.start[1],
                                  end[0], end[1])
        self.rotation = False
        return True

    def do_rotation(self, camera, window):
        if not self.rotation: return False
        new_end_screen = glfw.GetCursorPos(window)
        window_size = glfw.GetWindowSize(window)

        norm_x = -1.0 + 2.0 * new_end_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * new_end_screen[1] / window_size[1]
        new_end = (norm_x, norm_y)

        camera.update_orientation(self.start[0],
                                  self.start[1],
                                  new_end[0],
                                  new_end[1])
        self.start = new_end
        return True

