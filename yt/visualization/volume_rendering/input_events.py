# encoding: utf-8
"""
Input event handlers for Interactive Data Visualization

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# This is a part of the experimental Interactive Data Visualization

from collections import defaultdict, namedtuple
from yt.utilities.math_utils import \
    get_perspective_matrix, \
    get_orthographic_matrix
import OpenGL.GL as GL
import cyglfw3 as glfw
import numpy as np
import matplotlib.cm as cm
import random

event_registry = {}

GLFWEvent = namedtuple("GLFWEvent", ['window', 'key', 'scancode', 'action',
                       'mods', 'width', 'height'])

class EventCollection(object):
    '''Class handling mouse and keyboard events occurring in IDV
    
    Parameters
    ----------
    scene : :class:`yt.visualization.volume_rendering.interactive_vr.SceneGraph`
        A current scene object used in the IDV

    camera : :class:`yt.visualization.volume_rendering.interactive_vr.IDVCamera`
        A current camera object used in the IDV
    
    '''
    def __init__(self, scene, camera):
        self.key_callbacks = defaultdict(list)
        self.mouse_callbacks = defaultdict(list)
        self.framebuffer_callbacks = []
        self.render_events = []
        self.camera = camera
        self.scene = scene
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
    '''Close main window'''
    glfw.SetWindowShouldClose(event.window, True)

@register_event("zoomin")
def zoomin(event_coll, event):
    '''Zoom in the camera'''
    camera = event_coll.camera
    camera.position -= 0.05 * (camera.position - camera.focus) / \
                np.linalg.norm(camera.position - camera.focus)
    return True

@register_event("zoomout")
def zoomout(event_coll, event):
    '''Zoom out the camera'''
    camera = event_coll.camera
    camera.position += 0.05 * (camera.position - camera.focus) / \
        np.linalg.norm(camera.position - camera.focus)
    return True

@register_event("camera_orto")
def camera_orto(event_coll, event):
    '''Change camera to orthographic projection'''
    camera = event_coll.camera
    if camera.proj_func == get_orthographic_matrix:
        return False
    camera.proj_func = get_orthographic_matrix
    camera.fov = np.tan(np.radians(camera.fov) / 2.0)
    return True

@register_event("camera_proj")
def camera_proj(event_coll, event):
    '''Change camera to perspective projection'''
    camera = event_coll.camera
    if camera.proj_func == get_perspective_matrix:
        return False
    camera.proj_func = get_perspective_matrix
    camera.fov = np.degrees(np.arctan(camera.fov) * 2.0)
    return True


@register_event("shader_max")
def shader_max(event_coll, event):
    '''Use maximum intensity shader'''
    print("Changing shader to max(intensity)")
    scene = event_coll.scene
    for coll in scene.collections:
        coll.set_shader("default.v")
        coll.set_shader("max_intensity.f")
    scene.set_shader("passthrough.v")
    scene.set_shader("apply_colormap.f")
    for collection in scene.collections:
        collection.set_fields_log(True)
    scene.update_minmax()
    GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
    GL.glBlendEquation(GL.GL_MAX)
    return True

@register_event("shader_proj")
def shader_proj(event_coll, event):
    '''Use projection shader'''
    print("Changing shader to projection")
    scene = event_coll.scene
    for coll in scene.collections:
        coll.set_shader("default.v")
        coll.set_shader("projection.f")
    scene.set_shader("passthrough.v")
    scene.set_shader("apply_colormap.f")
    for collection in scene.collections:
        collection.set_fields_log(False)
    scene.update_minmax()
    GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
    GL.glBlendEquation(GL.GL_FUNC_ADD)
    return True

@register_event("shader_test")
def shader_test(event_coll, event):
    """Use transfer function shader"""
    print("Changing shader to projection")
    scene = event_coll.scene
    for coll in scene.collections:
        coll.set_shader("default.v")
        coll.set_shader("transfer_function.f")
    scene.set_shader("passthrough.v")
    scene.set_shader("noop.f")
    for collection in scene.collections:
        collection.set_fields_log(True)
    #scene.update_minmax()
    # https://www.opengl.org/sdk/docs/man/html/glBlendFunc.xhtml
    GL.glBlendEquationSeparate(GL.GL_FUNC_ADD, GL.GL_FUNC_ADD)
    GL.glBlendFuncSeparate(GL.GL_ONE, GL.GL_ONE_MINUS_SRC_ALPHA, GL.GL_ONE, GL.GL_ZERO)
    return True

@register_event("shader_lines")
def shader_lines(event_coll, event):
    print("Changing shader to projection")
    scene = event_coll.scene
    for coll in scene.collections:
        coll.set_shader("default.v")
        coll.set_shader("drawlines.f")
    scene.set_shader("passthrough.v")
    scene.set_shader("noop.f")
    for collection in scene.collections:
        collection.set_fields_log(True)
    #scene.update_minmax()
    # https://www.opengl.org/sdk/docs/man/html/glBlendFunc.xhtml
    #GL.glBlendFunc(GL.GL_ONE, GL.GL_DST_ALPHA)
    #GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
    #GL.glBlendFunc(GL.GL_ONE_MINUS_SRC_ALPHA, GL.GL_SRC_ALPHA)
    GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
    GL.glBlendEquation(GL.GL_MAX)
    return True


@register_event("cmap_cycle")
def cmap_cycle(event_coll, event):
    """Change colormap"""
    cmap = ['arbre', 'algae', 'kamae', 'viridis', 'inferno', 'magma']
    cmap = cm.get_cmap(random.choice(cmap))
    event_coll.camera.cmap = np.array(cmap(np.linspace(0, 1, 256)),
        dtype=np.float32)
    event_coll.camera.cmap_new = True
    print("Setting colormap to {}".format(cmap.name))
    return True

@register_event("cmap_max_up")
def cmap_max_up(event_coll, event):
    """Increase upper bound of colormap"""
    if event_coll.camera.cmap_log:
        event_coll.camera.cmap_max += 0.5
    else:
        event_coll.camera.cmap_max *= 2.0
    return True

@register_event("cmap_max_down")
def cmap_max_down(event_coll, event):
    """Decrease upper bound of colormap"""
    if event_coll.camera.cmap_log:
        event_coll.camera.cmap_max -= 0.5
    else:
        event_coll.camera.cmap_max *= 0.5
    return True

@register_event("cmap_min_up")
def cmap_min_up(event_coll, event):
    """Increase lower bound of colormap"""
    if event_coll.camera.cmap_log:
        event_coll.camera.cmap_min += 0.5
    else:
        event_coll.camera.cmap_min *= 2.0
    return True

@register_event("cmap_min_down")
def cmap_min_down(event_coll, event):
    """Decrease lower bound of colormap"""
    if event_coll.camera.cmap_log:
        event_coll.camera.cmap_min -= 0.5
    else:
        event_coll.camera.cmap_min *= 0.5
    return True

@register_event("cmap_toggle_log")
def cmap_toggle_log(event_coll, event):
    """Switch between linear and logarithmic scales"""
    if event_coll.scene.data_logged:
        print("Data is logged already, can't toggle scale to linear")
        return False

    if not event_coll.camera.cmap_log:
        event_coll.camera.cmap_max = np.log10(event_coll.camera.cmap_max)
        event_coll.camera.cmap_min = np.log10(event_coll.camera.cmap_min)
    else:
        event_coll.camera.cmap_max = 10.0 ** event_coll.camera.cmap_max
        event_coll.camera.cmap_min = 10.0 ** event_coll.camera.cmap_min
    event_coll.camera.cmap_log = not event_coll.camera.cmap_log
    return True

@register_event("closeup")
def closeup(event_coll, event):
    """Change camera position to (0.01, 0.01, 0.01)"""
    event_coll.camera.position = (0.01, 0.01, 0.01)
    return True

@register_event("reset")
def reset(event_coll, event):
    """Change camera position to (-1.0, -1.0, -1.0)"""
    event_coll.camera.position = (-1.0, -1.0, -1.0)
    return True

@register_event("print_limits")
def print_limits(event_coll, event):
    """Print debug info about scene and camera"""
    print(event_coll.scene.min_val, event_coll.scene.max_val)
    print(event_coll.camera.cmap_min, event_coll.camera.cmap_max,
          event_coll.camera.cmap_log)
    return False

@register_event("debug_buffer")
def debug_buffer(event_coll, event):
    """Print debug info about framebuffer"""
    buffer = event_coll.scene._retrieve_framebuffer()
    print(buffer.min(), buffer.max())

@register_event("print_help")
def print_help(event_coll, event):
    """Print this help"""
    key_map = {}
    for key in (a for a in dir(glfw) if a.startswith("KEY")):
        key_map[glfw.__dict__.get(key)] = key[4:]
    for cb in (f for f in sorted(event_coll.key_callbacks)
               if isinstance(f, tuple)):
        print("%s - %s" % (key_map[cb[0]],
                           event_coll.key_callbacks[cb][0].__doc__))
    return False

@register_event("nplane_closer")
def nplane_closer(event_coll, event):
    print("nearplane", event_coll.camera.near_plane)
    event_coll.camera.near_plane /= 2.0
    return True

@register_event("nplane_further")
def nplane_further(event_coll, event):
    print("nearplane", event_coll.camera.near_plane)
    event_coll.camera.near_plane *= 2.0
    return True

class MouseRotation(object):
    '''Class translating mouse movements to positions in OpenGL scene's coordinates'''
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

class BlendFuncs(object):
    '''Class allowing to switch between different GL blending functions'''
    possibilities = (
        "GL_ZERO", "GL_ONE", "GL_SRC_COLOR", "GL_ONE_MINUS_SRC_COLOR",
        "GL_DST_COLOR", "GL_ONE_MINUS_DST_COLOR", "GL_SRC_ALPHA",
        "GL_ONE_MINUS_SRC_ALPHA", "GL_DST_ALPHA", "GL_ONE_MINUS_DST_ALPHA",
        "GL_CONSTANT_COLOR", "GL_ONE_MINUS_CONSTANT_COLOR",
        "GL_CONSTANT_ALPHA", "GL_ONE_MINUS_CONSTANT_ALPHA")
    source_i = 0
    dest_i = 0

    def next_source(self, event_coll, event):
        self.source_i = (self.source_i + 1) % len(self.possibilities)
        s = getattr(GL, self.possibilities[self.source_i])
        d = getattr(GL, self.possibilities[self.dest_i])
        print("Setting source to %s and dest to %s" %
              (self.possibilities[self.source_i], 
               self.possibilities[self.dest_i]))
        GL.glBlendFunc(s, d)
        return True

    def next_dest(self, event_coll, event):
        self.dest_i = (self.dest_i + 1) % len(self.possibilities)
        s = getattr(GL, self.possibilities[self.source_i])
        d = getattr(GL, self.possibilities[self.dest_i])
        print("Setting source to %s and dest to %s" %
              (self.possibilities[self.source_i],
               self.possibilities[self.dest_i]))
        GL.glBlendFunc(s, d)
        return True
