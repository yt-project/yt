import cyglfw3
import imgui
import contextlib

from imgui.integrations import compute_fb_scale
from imgui.opengl import ProgrammablePipelineRenderer

class Cyglfw3Renderer(ProgrammablePipelineRenderer):
    def __init__(self, window, attach_callbacks=True):
        super(Cyglfw3Renderer, self).__init__()
        self.window = window

        io = imgui.get_io()
        io.display_size = cyglfw3.GetFramebufferSize(self.window)

        self._map_keys()
        self._gui_time = None

    def _attach_callbacks(self, callbacks):
        # Note that this *un-attaches* the previous attached items.
        cyglfw3.SetKeyCallback(self.window, self.keyboard_callback)
        cyglfw3.SetCursorPosCallback(self.window, self.cursor_callback)
        cyglfw3.SetMouseButtonCallback(self.window, self.mouse_callback)
        cyglfw3.SetWindowSizeCallback(self.window, self.resize_callback)
        cyglfw3.SetCharCallback(self.window, self.char_callback)
        cyglfw3.SetScrollCallback(self.window, self.scroll_callback)

        self.callbacks = callbacks

    def _map_keys(self):
        key_map = self.io.key_map

        key_map[imgui.KEY_TAB] = cyglfw3.KEY_TAB
        key_map[imgui.KEY_LEFT_ARROW] = cyglfw3.KEY_LEFT
        key_map[imgui.KEY_RIGHT_ARROW] = cyglfw3.KEY_RIGHT
        key_map[imgui.KEY_UP_ARROW] = cyglfw3.KEY_UP
        key_map[imgui.KEY_DOWN_ARROW] = cyglfw3.KEY_DOWN
        key_map[imgui.KEY_PAGE_UP] = cyglfw3.KEY_PAGE_UP
        key_map[imgui.KEY_PAGE_DOWN] = cyglfw3.KEY_PAGE_DOWN
        key_map[imgui.KEY_HOME] = cyglfw3.KEY_HOME
        key_map[imgui.KEY_END] = cyglfw3.KEY_END
        key_map[imgui.KEY_DELETE] = cyglfw3.KEY_DELETE
        key_map[imgui.KEY_BACKSPACE] = cyglfw3.KEY_BACKSPACE
        key_map[imgui.KEY_ENTER] = cyglfw3.KEY_ENTER
        key_map[imgui.KEY_ESCAPE] = cyglfw3.KEY_ESCAPE
        key_map[imgui.KEY_A] = cyglfw3.KEY_A
        key_map[imgui.KEY_C] = cyglfw3.KEY_C
        key_map[imgui.KEY_V] = cyglfw3.KEY_V
        key_map[imgui.KEY_X] = cyglfw3.KEY_X
        key_map[imgui.KEY_Y] = cyglfw3.KEY_Y
        key_map[imgui.KEY_Z] = cyglfw3.KEY_Z

    def keyboard_callback(self, window, key, scancode, action, mods):
        # perf: local for faster access
        io = imgui.get_io()

        if action == cyglfw3.PRESS:
            io.keys_down[key] = True
        elif action == cyglfw3.RELEASE:
            io.keys_down[key] = False

        io.key_ctrl = (
            io.keys_down[cyglfw3.KEY_LEFT_CONTROL] or
            io.keys_down[cyglfw3.KEY_RIGHT_CONTROL]
        )

        io.key_alt = (
            io.keys_down[cyglfw3.KEY_LEFT_ALT] or
            io.keys_down[cyglfw3.KEY_RIGHT_ALT]
        )

        io.key_shift = (
            io.keys_down[cyglfw3.KEY_LEFT_SHIFT] or
            io.keys_down[cyglfw3.KEY_RIGHT_SHIFT]
        )

    def char_callback(self, window, char):
        io = imgui.get_io()

        if 0 < char < 0x10000:
            io.add_input_character(char)

    def resize_callback(self, window, width, height):
        self.io.display_size = width, height

    def cursor_callback(self, *args, **kwargs):
        pass

    def mouse_callback(self, *args, **kwargs):
        io = imgui.get_io()

    def scroll_callback(self, window, x_offset, y_offset):
        io = imgui.get_io()
        io.mouse_wheel = y_offset

    def process_inputs(self):
        io = imgui.get_io()

        window_size = cyglfw3.GetWindowSize(self.window)
        fb_size = cyglfw3.GetFramebufferSize(self.window)

        io.display_size = window_size
        io.display_fb_scale = compute_fb_scale(window_size, fb_size)
        io.delta_time = 1.0/60

        if cyglfw3.GetWindowAttrib(self.window, cyglfw3.FOCUSED):
            io.mouse_pos = cyglfw3.GetCursorPos(self.window)
        else:
            io.mouse_pos = -1, -1

        io.mouse_down[0] = cyglfw3.GetMouseButton(self.window, 0)
        io.mouse_down[1] = cyglfw3.GetMouseButton(self.window, 1)
        io.mouse_down[2] = cyglfw3.GetMouseButton(self.window, 2)

        current_time = cyglfw3.GetTime()

        if self._gui_time:
            io.delta_time = current_time - self._gui_time
        else:
            io.delta_time = 1. / 60.

        self._gui_time = current_time

class SimpleGUI:

    renderer = None
    callbacks = None
    context = None
    window = None
    draw = False

    def __init__(self, window, callbacks):
        self.window = window
        self.callbacks = callbacks
        self.context = imgui.create_context()
        self.renderer = Cyglfw3Renderer(self.window, attach_callbacks = False)

        # Setup the IO
        io = imgui.get_io()
        io.fonts.add_font_default()
        io.display_size = cyglfw3.glfw3.GetWindowSize(self.window)

    def init_frame(self):
        self.renderer.process_inputs()
        imgui.new_frame()

    def finalize_frame(self):
        imgui.render()

    def draw_frame(self):
        self.renderer.render(imgui.get_draw_data())

    def run(self, scene, callbacks):
        if imgui.begin_main_menu_bar():
            if imgui.begin_menu("File", True):
                clicked_quit, selected_quit = imgui.menu_item(
                    "Quit", "Cmd+Q", False, True
                )
                imgui.end_menu()
        imgui.end_main_menu_bar()
        imgui.show_test_window()
        imgui.begin("Custom window", True)
        imgui.text("Bar")
        imgui.text_colored("Eggs", 0.2, 1.0, 0.0)
        imgui.end()
