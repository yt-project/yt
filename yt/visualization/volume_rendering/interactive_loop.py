import cyglfw3 as glfw
import numpy as np
from OpenGL.GL import glViewport

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

def key_callback(window, key, scancode, action, mods):
    global draw, c
    if key == glfw.KEY_ESCAPE and action == glfw.PRESS:
        glfw.SetWindowShouldClose(window, True)
    if key == glfw.KEY_W and action == glfw.PRESS:
        c.position -= (c.position - c.focus) / np.linalg.norm(c.position - c.focus)
        draw = True
    if key == glfw.KEY_S and action == glfw.PRESS:
        c.position += (c.position - c.focus) / np.linalg.norm(c.position - c.focus)
        draw = True

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
    glfw.SetKeyCallback(window, key_callback)
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

        if draw:
            scene.set_camera(c)
            scene.render()
            draw = False
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
