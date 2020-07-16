import imgui
from imgui.integrations.pyglet import PygletRenderer
import contextlib

class SimpleGUI:

    renderer = None
    callbacks = None
    context = None
    window = None
    draw = False

    def __init__(self, window):
        self.window = window
        self.context = imgui.create_context()
        self.renderer = PygletRenderer(window)

    def render(self, scene):
        imgui.new_frame()
        changed = False
        if scene is not None:
            for i, component in enumerate(scene.components):
                imgui.begin("Component {}".format(i))
                _, component.cmap_log = imgui.checkbox("Take log", component.cmap_log)
                changed = changed or _
                imgui.end()
            self.window._do_update = self.window._do_update or changed
        imgui.render()
        self.renderer.render(imgui.get_draw_data())

    @property
    def mouse_event_handled(self):
        return self.renderer.io.want_capture_mouse

    @property
    def keyboard_event_handled(self):
        return self.renderer.io.want_capture_keyboard
