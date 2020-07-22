import array
import ctypes
import itertools
import pyglet
import numpy as np
from math import ceil, floor
import matplotlib.pyplot as plt

import imgui
from imgui.integrations.pyglet import create_renderer
import contextlib
from ..image_writer import write_bitmap

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
        self.snapshot_count = 0
        self.snapshot_format = r"snap_{count:04d}.png"

    def render(self, scene):
        imgui.new_frame()
        changed = False
        if scene is not None:
            imgui.style_colors_classic()
            imgui.begin("Scene")
            imgui.text("Filename Template:")
            _, self.snapshot_format = imgui.input_text("", self.snapshot_format, 256)
            if imgui.button("Save Snapshot"):
                # Call render again, since we're in the middle of overlaying
                # some stuff and we want a clean scene snapshot
                scene.render()
                write_bitmap(scene.image[:,:,:3], self.snapshot_format.format(count = self.snapshot_count))
                self.snapshot_count += 1
            imgui.text("Camera Position")
            pos = scene.camera.position
            _, values = imgui.input_float3("", pos[0], pos[1], pos[2], flags =
                                           imgui.INPUT_TEXT_ENTER_RETURNS_TRUE)
            if _:
                scene.camera.position = np.array(values)
                scene.camera._update_matrices()
            changed = changed or _
            #imgui.show_style_editor()
            for i, element in enumerate(scene):
                if imgui.tree_node("element {}: {}".format(i + 1, element.name)):
                    changed = changed or element.render_gui(imgui, self.renderer)
                    imgui.tree_pop()
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

# https://bitbucket.org/pyglet/pyglet/issues/143/support-numpy-in-more-places
def get_numpy_data(arr):
     """
     :param arr: numpy array of float32
     :return: ctypes array of float32
     """
     # Accept any contiguous array of float32
     assert arr.flags["C_CONTIGUOUS"] or arr.flags["F_CONTIGUOUS"]
     assert arr.dtype == np.uint8
     return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8 * arr.size))[0]

class WindowWidget:
    def __init__(self, window, x0, y0, width, height, draw_border = False):
        self.draw_border = draw_border
        self.window = window
        self.x0 = x0
        self.y0 = y0
        self.width = width
        self.height = height
        self.x1 = x0 + width
        self.y1 = y0 + height
        self.setup()

    def local_pos(self, x, y):
        return x - self.x0, y - self.y0

    def local_scale(self, x, y):
        return (x - self.x0)/self.width, (y - self.y0)/self.height

    def do_draw(self):
        pass

    def on_draw(self):
        self.do_draw()
        if self.draw_border:
            pyglet.graphics.draw(4, pyglet.gl.GL_LINE_LOOP,
                ('vertices2f', [self.x0, self.y0,
                         self.x0 + self.width, self.y0,
                         self.x0 + self.width, self.y0 + self.height,
                         self.x0, self.y0 + self.height]),
                ('colors3bn', [255 for _ in range(4 * 3)])
            )

class TransferFunctionWidget(WindowWidget):

    def setup(self):
        self.active = 'red'
        self.vals = {}
        self.colors = {}
        self.draw = {}
        self.N_bins = 256

        viridis = np.array(plt.get_cmap("viridis").colors)
        for i, color in enumerate(('red', 'green', 'blue')):
            self.vals[color] = np.vstack([np.mgrid[0.0:1.0:1j*self.N_bins],
                                     np.mgrid[0.0:1.0:self.N_bins*1j]]).copy(order="F")
            self.vals[color][1,:] = viridis[:,i]
            self.colors[color] = np.zeros((self.N_bins, 3), dtype="u1").copy(order="C")
            self.colors[color][:,i] = 255
            self.update_draw(color)
        self.vals['alpha'] = self.vals['red'].copy()
        self.vals['alpha'][1,:] = 1.0
        self.colors['alpha'] = np.zeros((self.N_bins, 3), dtype="u1").copy(order="C") + 255
        self.update_draw('alpha')

    def update_draw(self, color):
        self.draw[color] = self.vals[color] \
            * np.array([self.width, self.height])[:, None] \
            + np.array([self.x0, self.y0])[:, None]


    def compute_bin_range(self, x, dx):
        left_bound = floor(min(x + dx, x) * self.N_bins / self.width)
        right_bound = ceil(max(x + dx, x) * self.N_bins / self.width)
        return (max(min(_, self.N_bins -1), 0) for _ in (left_bound, right_bound))

    def compute_val_range(self, y, dy):
        yv1 = (y / self.height)
        yv2 = (y + dy) / self.height
        return (max(min(_, 1.0), 0.0) for _ in (yv1, yv2))

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if x < self.x0 or x > self.x1 or y < self.y0 or y > self.y1:
            return False
        xb1, xb2 = self.compute_bin_range(x - self.x0, dx)
        yv1, yv2 = self.compute_val_range(y - self.y0, dy)
        if buttons & pyglet.window.mouse.RIGHT:
            if modifiers & pyglet.window.key.MOD_SHIFT:
                yv1 = yv2 = 0.0
            else:
                yv1 = yv2 = 1.0

        if dx < 0: yv2, yv1 = yv1, yv2
        self.vals[self.active][1,xb1:xb2] = np.mgrid[yv1:yv2:(xb2 - xb1)*1j]
        self.update_draw(self.active)
        return True

    def on_key_press(self, symbol, modifiers):
        if symbol == pyglet.window.key.R:
            self.active = "red"
            return True
        elif symbol == pyglet.window.key.G:
            self.active = "green"
            return True
        elif symbol == pyglet.window.key.B:
            self.active = "blue"
            return True
        elif symbol == pyglet.window.key.A:
            self.active = "alpha"
            return True
        return False

    def do_draw(self):
        for color in ('red', 'green', 'blue', 'alpha'):
            pyglet.graphics.draw(self.N_bins, pyglet.gl.GL_LINE_STRIP,
                ('vertices2f', self.draw[color].ravel(order="F")),
                ('colors3Bn', self.colors[color].ravel(order="C")))

class TransferFunctionImage(WindowWidget):
    def __init__(self, window, x0, y0, width, height, tf_widget, draw_border = False):
        self.tf_widget = tf_widget
        super(TransferFunctionImage, self).__init__(window, x0, y0, width, height, draw_border = draw_border)

    def create_image_data(self):
        ii = np.vstack([(self.tf_widget.vals[c][1,:] * 255).astype("u1")
                         for c in ('red', 'green', 'blue', 'alpha')]).copy(order="F")
        self.image_data = pyglet.image.ImageData(
            self.tf_widget.N_bins, 1,
            'RGBA', get_numpy_data(ii),
            3*self.tf_widget.N_bins
        )

    def setup(self):
        self.create_image_data()
        self.sprite = pyglet.sprite.Sprite(self.image_data, self.x0, self.y0)
        self.sprite.update(scale_x = self.width / self.tf_widget.N_bins,
                           scale_y = self.height)

    def on_key_press(self, *args, **kwargs):
        pass

    def on_mouse_drag(self, *args, **kwargs):
        pass

    def update_sprite(self):
        self.create_image_data()
        self.sprite.image = self.image_data

    def do_draw(self):
        self.update_sprite()
        self.sprite.draw()

class TransferFunctionWindow(pyglet.window.Window):
    def __init__(self, *args, **kwargs):
        super(TransferFunctionWindow, self).__init__(*args, **kwargs)
        self.widgets = []
        print("original", self.width, self.height)

    def on_draw(self):
        self.switch_to()
        self.clear()
        for w in self.widgets:
            w.on_draw()
        self.flip()

    def on_resize(self, *args, **kwargs):
        print("Resize", args, kwargs)
        super(TransferFunctionWindow, self).on_resize(*args, **kwargs)

    def on_key_press(self, symbol, modifiers):
        for w in self.widgets:
            if w.on_key_press(symbol, modifiers):
                break
        else:
            super(TransferFunctionWindow, self).on_key_press(symbol, modifiers)

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        for w in self.widgets:
            if w.on_mouse_drag(x, y, dx, dy, buttons, modifiers):
                break
