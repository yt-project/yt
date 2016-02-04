import cyglfw3 as glfw
import numpy as np
import OpenGL.GL as GL
from collections import defaultdict

from interactive_vr import BlockCollection, SceneGraph, TrackballCamera


class Events(object):
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
            draw = f(window, width, height) or draw
        self.draw = self.draw or draw

    def __call__(self, window):
        draw = False
        for f in self.render_events:
            draw = f(self.camera, window) or draw
        self.draw = self.draw or draw

    def add_default(self, func):
        self.render_events.append(func)

    def add_key_callback(self, func, key, action = "press", mods = None):
        self._add_callback(self.key_callbacks, func, key, action, mods)

    def add_mouse_callback(self, func, key, action = "press", mods = None):
        self._add_callback(self.mouse_callbacks, func, key, action, mods)

    def add_framebuffer_callback(self, func):
        self.framebuffer_callbacks.append(func)

    def _add_callback(self, d, func, key, action, mods):
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

def framebuffer_size_callback(window, width, height):
    GL.glViewport(0, 0, width, height)
    return True

def close_window(camera, window, key, scancode, action, mods):
    glfw.SetWindowShouldClose(window, True)

def zoomin(camera, window, key, scancode, action, mods):
    camera.position -= 0.05 * (camera.position - camera.focus) / \
                np.linalg.norm(camera.position - camera.focus)
    print camera.position, camera.focus
    return True

def zoomout(camera, window, key, scancode, action, mods):
    camera.position += 0.05 * (camera.position - camera.focus) / \
        np.linalg.norm(camera.position - camera.focus)
    print camera.position, camera.focus
    return True

def closeup(camera, window, key, scancode, action, mods):
    camera.position = (0.01, 0.01, 0.01)
    return True

def reset(camera, window, key, scancode, action, mods):
    camera.position = (-1.0, -1.0, -1.0)
    return True

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
        
class RenderingContext(object):
    def __init__(self, width = 800, height = 600, title = "vol_render"):
        glfw.Init()
        glfw.WindowHint(glfw.CONTEXT_VERSION_MAJOR, 3)
        glfw.WindowHint(glfw.CONTEXT_VERSION_MINOR, 3)
        glfw.WindowHint(glfw.OPENGL_FORWARD_COMPAT, True)
        glfw.WindowHint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
        self.window = glfw.CreateWindow(width, height, title)
        if not self.window:
            glfw.Terminate()
            exit()

        glfw.MakeContextCurrent(self.window)
        GL.glClearColor(0.0, 0.0, 0.0, 0.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        glfw.SwapBuffers(self.window)
        glfw.PollEvents()

    def start_loop(self, scene, camera):
        scene.set_camera(camera)
        scene.add_shader_from_file("max_intensity_frag.glsl")
        camera.compute_matrices()
        frame_start = glfw.GetTime()
        fps_start = glfw.GetTime()
        f = 0
        N = 10.0
        print "Starting rendering..."
        callbacks = Events(camera)
        callbacks.add_key_callback(close_window, "escape")
        callbacks.add_key_callback(zoomin, "w")
        callbacks.add_key_callback(zoomout, "s")
        callbacks.add_key_callback(closeup, "z")
        callbacks.add_key_callback(reset, "r")
        mouse_callbacks = MouseRotation()
        callbacks.add_mouse_callback(mouse_callbacks.start_rotation,
            glfw.MOUSE_BUTTON_LEFT)
        callbacks.add_mouse_callback(mouse_callbacks.stop_rotation,
            glfw.MOUSE_BUTTON_LEFT, action="release")
        callbacks.add_framebuffer_callback(framebuffer_size_callback)
        callbacks.add_default(mouse_callbacks.do_rotation)
        glfw.SetFramebufferSizeCallback(self.window,
            callbacks.framebuffer_call)
        glfw.SetKeyCallback(self.window, callbacks.key_call)
        glfw.SetMouseButtonCallback(self.window, callbacks.mouse_call)
        callbacks.draw = True
        while not glfw.WindowShouldClose(self.window):
            callbacks(self.window)
            if callbacks.draw:
                camera.compute_matrices()
                scene.set_camera(camera)
                scene.render()
                callbacks.draw = False
                glfw.SwapBuffers(self.window)
            glfw.PollEvents()

        print "Finished rendering"
        glfw.Terminate()
