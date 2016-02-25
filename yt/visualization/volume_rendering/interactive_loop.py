import cyglfw3 as glfw
import numpy as np
import OpenGL.GL as GL
from .input_events import EventCollection, MouseRotation

from interactive_vr import BlockCollection, SceneGraph, TrackballCamera

class RenderingContext(object):
    should_quit = False
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

    def setup_loop(self, scene, camera):
        scene.set_camera(camera)
        scene.update_minmax()
        camera.compute_matrices()
        frame_start = glfw.GetTime()
        fps_start = glfw.GetTime()
        print "Starting rendering..."
        callbacks = EventCollection(scene, camera)
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
        callbacks.add_key_callback("print_limits", "h")
        callbacks.add_key_callback("debug_buffer", "d")
        callbacks.add_key_callback("cmap_max_up", "right_bracket")
        callbacks.add_key_callback("cmap_max_down", "left_bracket")
        callbacks.add_key_callback("cmap_min_up", "semicolon")
        callbacks.add_key_callback("cmap_min_down", "apostrophe")
        callbacks.add_key_callback("cmap_toggle_log", "l")
        mouse_callbacks = MouseRotation()
        callbacks.add_mouse_callback(mouse_callbacks.start_rotation,
            glfw.MOUSE_BUTTON_LEFT)
        callbacks.add_mouse_callback(mouse_callbacks.stop_rotation,
            glfw.MOUSE_BUTTON_LEFT, action="release")
        callbacks.add_framebuffer_callback("framebuffer_size_callback")
        callbacks.add_render_callback(mouse_callbacks.do_rotation)
        glfw.SetFramebufferSizeCallback(self.window,
            callbacks.framebuffer_call)
        glfw.SetKeyCallback(self.window, callbacks.key_call)
        glfw.SetMouseButtonCallback(self.window, callbacks.mouse_call)
        callbacks.draw = True
        return callbacks

    def start_loop(self, scene, camera):
        callbacks = self.setup_loop(scene, camera)
        for i in self(scene, camera, callbacks):
            pass

    def __call__(self, scene, camera, callbacks):
        while not glfw.WindowShouldClose(self.window) or self.should_quit:
            callbacks(self.window)
            if callbacks.draw:
                camera.compute_matrices()
                scene.set_camera(camera)
                scene.render()
                glfw.SwapBuffers(self.window)
                callbacks.draw = False
            glfw.PollEvents()
            yield self
        glfw.Terminate()
