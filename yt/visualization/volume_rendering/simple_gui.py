import imgui
from imgui.integrations.pyglet import create_renderer
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
        self.renderer = create_renderer(window)
        self.expanded = {}

    def render(self, scene):
        imgui.new_frame()
        changed = False
        if scene is not None:
            imgui.begin("Scene")
            elements = scene.components + scene.annotations
            for i, element in enumerate(sorted(elements, key = lambda a: a.priority)):
                self.expanded.setdefault(id(element), True)
                expanded, visible = imgui.collapsing_header("element {}: {}".format(i + 1, element.name))
                self.expanded[id(element)] = expanded
                if self.expanded[id(element)]:
                    changed = changed or element.render_gui(imgui, self.renderer)
            self.window._do_update = self.window._do_update or changed
            imgui.end()
        imgui.render()
        self.renderer.render(imgui.get_draw_data())

    @property
    def mouse_event_handled(self):
        return self.renderer.io.want_capture_mouse

    @property
    def keyboard_event_handled(self):
        return self.renderer.io.want_capture_keyboard
