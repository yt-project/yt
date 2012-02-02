from yt.mods import *
import numpy as na

from yt.utilities.answer_testing.output_tests import \
    YTStaticOutputTest, RegressionTestException
from yt.funcs import ensure_list


class VolumeRenderingInconsistent(RegressionTestException):
    pass


class VolumeRenderingConsistency(YTStaticOutputTest):
    name = "volume_rendering_consistency"

    def run(self):
        c = (self.pf.domain_right_edge + self.pf.domain_left_edge) / 2.
        W = na.sqrt(3.) * (self.pf.domain_right_edge - \
            self.pf.domain_left_edge)
        N = 512
        n_contours = 5
        cmap = 'algae'
        field = 'Density'
        mi, ma = self.pf.h.all_data().quantities['Extrema'](field)[0]
        mi, ma = na.log10(mi), na.log10(ma)
        contour_width = (ma - mi) / 100.
        L = na.array([1.] * 3)
        tf = ColorTransferFunction((mi - 2, ma + 2))
        tf.add_layers(n_contours, w=contour_width,
                      col_bounds=(mi * 1.001, ma * 0.999),
                      colormap=cmap, alpha=na.logspace(-1, 0, n_contours))
        cam = self.pf.h.camera(c, L, W, (N, N), transfer_function=tf,
            no_ghost=True)
        image = cam.snapshot()
        # image = cam.snapshot('test_rendering_%s.png'%field)
        self.result = image

    def compare(self, old_result):
        # Compare the deltas; give a leeway of 1e-8
        delta = na.nanmax(na.abs(self.result - old_result) /
                                 (self.result + old_result))
        if delta > 1e-9: raise VolumeRenderingInconsistent()
