
from imgui.integrations.cyglfw3 import Cyglfw3Renderer
import cyglfw3
import imgui
import contextlib

class SimpleGUI:

    renderer = None
    window = None
    io = None
    _initialized = False
    draw = False

    def __init__(self, window):
        self.window = window

    def setup_renderer(self):
        if self.renderer is not None: return
        self.context = imgui.create_context()
        self.renderer = Cyglfw3Renderer(self.window, attach_callbacks = False)
        self.io = imgui.get_io()
        self.io.fonts.add_font_default()
        self.io.display_size = cyglfw3.glfw3.GetWindowSize(self.window)

    def init_frame(self):
        cyglfw3.PollEvents()
        imgui.new_frame()
        self.renderer.process_inputs()

    def finalize_frame(self):
        imgui.render()
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

    def respond_to_inputs(self, callbacks):
        if not self._initialized:
            # this is where we update callbacks to have ours go first.
            def _mouse_input(callbacks):
                def check_mouse_input(*args, **kwargs):
                    print("Checking mouse", self.io.want_capture_mouse)
                    if not self.io.want_capture_mouse:
                        callbacks.mouse_call(*args, **kwargs)
                    else:
                        self.draw = True
                return check_mouse_input
            def _keyboard_input(callbacks):
                def check_keyboard_input(*args, **kwargs):
                    print("Checking keyboard", self.io.want_capture_keyboard)
                    if not self.io.want_capture_keyboard:
                        callbacks.key_call(*args, **kwargs)
                    else:
                        self.draw = True
                return check_keyboard_input
            callbacks.framebuffer_callbacks.insert(0, self.renderer.resize_callback)
            # Re-insert the ones we need for the renderer
            cyglfw3.SetMouseButtonCallback(self.window, _mouse_input(callbacks))
            cyglfw3.SetCharCallback(self.window, self.renderer.char_callback)
            cyglfw3.SetScrollCallback(self.window, self.renderer.scroll_callback)
            cyglfw3.SetCursorPosCallback(self.window, self.renderer.cursor_callback)
            self._initialized = True

    def process_keyboard_inputs(self):
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

        #print("io.mouse_down: {} {} {}".format(io.mouse_down[0], io.mouse_down[1], io.mouse_down[2]))
        io.mouse_down[0] = cyglfw3.GetMouseButton(self.window, 0)
        io.mouse_down[1] = cyglfw3.GetMouseButton(self.window, 1)
        io.mouse_down[2] = cyglfw3.GetMouseButton(self.window, 2)

        current_time = cyglfw3.GetTime()

        if self._gui_time:
            self.io.delta_time = current_time - self._gui_time
        else:
            self.io.delta_time = 1. / 60.

        self._gui_time = current_time
