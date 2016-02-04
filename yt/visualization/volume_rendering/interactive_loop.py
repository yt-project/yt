import cyglfw3 as glfw
import numpy as np
from OpenGL.GL import glViewport
from collections import defaultdict

from interactive_vr import BlockCollection, SceneGraph, Camera

draw = True
start = None
rotation = False
window = None

c = None

def mouse_button_callback(window, button, action, mods):
    global draw, start, rotation, c

    if (button == glfw.MOUSE_BUTTON_LEFT and action == glfw.PRESS):
        start_screen = glfw.GetCursorPos(window) # Screen coordinates
        window_size = glfw.GetWindowSize(window)

        norm_x = -1.0 + 2.0 * start_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * start_screen[1] / window_size[1]
        start = (norm_x, norm_y)
        rotation = True

    if (button == glfw.MOUSE_BUTTON_LEFT and action == glfw.RELEASE):
        end_screen = glfw.GetCursorPos(window)
        window_size = glfw.GetWindowSize(window)

        norm_x = -1.0 + 2.0 * end_screen[0] / window_size[0]
        norm_y = 1.0 - 2.0 * end_screen[1] / window_size[1]
        end = (norm_x, norm_y)

        c.update_position( (end[0] - start[0]), (end[1] - start[1]))

        rotation = False
        draw = True

def framebuffer_size_callback(window, width, height):
    global draw
    glViewport(0, 0, width, height)
    draw = True

class KeyCallbacks(object):
    draw = True
    def __init__(self, camera):
        self.callbacks = defaultdict(list)
        self.camera = camera

    def __call__(self, window, key, scancode, action, mods):
        draw = False
        for f in self.callbacks[None] + self.callbacks[key, action, mods]:
            draw = f(self.camera, window, key, scancode, action, mods) or draw
        self.draw = draw

    def add_default(self, func):
        self.callbacks[None].append(func)

    def add(self, func, key, action = "press", mods = None):
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
        self.callbacks[key, action, mod].append(func)

def close_window(camera, window, key, scancode, action, mods):
    glfw.SetWindowShouldClose(window, True)

def zoomin(camera, window, key, scancode, action, mods):
    camera.position -= (camera.position - camera.focus) / \
                np.linalg.norm(camera.position - camera.focus)
    return True

def zoomout(camera, window, key, scancode, action, mods):
    camera.position += (camera.position - camera.focus) / \
        np.linalg.norm(camera.position - camera.focus)
    return True

def start_context():
    global window
    glfw.Init()

    glfw.WindowHint(glfw.CONTEXT_VERSION_MAJOR, 3)
    glfw.WindowHint(glfw.CONTEXT_VERSION_MINOR, 3)
    glfw.WindowHint(glfw.OPENGL_FORWARD_COMPAT, True)
    glfw.WindowHint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
    window = glfw.CreateWindow(800, 600, 'vol_render')
    if not window:
        glfw.Terminate()
        exit()

    glfw.MakeContextCurrent(window)
    glfw.SetFramebufferSizeCallback(window, framebuffer_size_callback)
    glfw.SetMouseButtonCallback(window, mouse_button_callback)

def start_loop(scene, camera):
    global draw, c, start
    c = camera
    scene.set_camera(camera)
    scene.add_shader_from_file("max_intensity_frag.glsl")
    frame_start = glfw.GetTime()
    fps_start = glfw.GetTime()
    f = 0
    N = 10.0
    print "Starting rendering..."
    callbacks = KeyCallbacks(camera)
    callbacks.add(close_window, "escape")
    callbacks.add(zoomin, "w")
    callbacks.add(zoomout, "s")
    glfw.SetKeyCallback(window, callbacks)

    while not glfw.WindowShouldClose(window):
        if rotation:
            new_end_screen = glfw.GetCursorPos(window)
            window_size = glfw.GetWindowSize(window)

            norm_x = -1.0 + 2.0 * new_end_screen[0] / window_size[0]
            norm_y = 1.0 - 2.0 * new_end_screen[1] / window_size[1]
            new_end = (norm_x, norm_y)

            c.update_position( (new_end[0] - start[0]),
                    (new_end[1] - start[1]))
            start = new_end
            draw = True

        if draw or callbacks.draw:
            scene.set_camera(c)
            scene.render()
            draw = callbacks.draw = False
        glfw.SwapBuffers(window)
        glfw.PollEvents()

        frame_start = glfw.GetTime()
        f += 1
        if f == N:
            print "FPS:", N / float(frame_start - fps_start)
            fps_start = glfw.GetTime()
            f = 0

    print "Finished rendering"
    glfw.Terminate()
