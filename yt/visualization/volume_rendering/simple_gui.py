import imgui
from imgui.integrations.pyglet import create_renderer
import contextlib

cmaps = ["arbre", "viridis", "magma", "doom"]

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

    def render(self, scene):
        imgui.new_frame()
        changed = False
        if scene is not None:
            for i, component in enumerate(scene.components):
                component.render_gui(i, imgui, self.renderer)
                continue
                imgui.begin("Component {}".format(i))
                _, component.cmap_log = imgui.checkbox("Take log", component.cmap_log)
                changed = changed or _
                _, cmap_index = imgui.listbox(
                    "Colormap", cmaps.index(component.colormap.colormap_name), cmaps)
                if _: component.colormap.colormap_name = cmaps[cmap_index]
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
