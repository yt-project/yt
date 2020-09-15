from functools import wraps

import numpy as np

filter_registry = {}


def apply_filter(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        args[0]._filters.append((f.__name__, (args, kwargs)))
        return args[0]

    return newfunc


class FixedResolutionBufferFilter:

    """
    This object allows to apply data transformation directly to
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    """

    def __init_subclass__(cls, *args, **kwargs):
        super().__init_subclass__(*args, **kwargs)
        filter_registry[cls.__name__] = cls

    def __init__(self, *args, **kwargs):
        pass

    def apply(self, buff):
        pass


class FixedResolutionBufferGaussBeamFilter(FixedResolutionBufferFilter):

    """
    This filter convolves
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer` with
    2d gaussian that is 'nbeam' pixels wide and has standard deviation
    'sigma'.
    """

    _filter_name = "gauss_beam"

    def __init__(self, nbeam=30, sigma=2.0):
        self.nbeam = nbeam
        self.sigma = sigma

    def apply(self, buff):
        from yt.utilities.on_demand_imports import _scipy

        hnbeam = self.nbeam // 2
        sigma = self.sigma

        l = np.linspace(-hnbeam, hnbeam, num=self.nbeam + 1)
        x, y = np.meshgrid(l, l)
        g2d = (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(
            -((x / sigma) ** 2 + (y / sigma) ** 2) / (2 * sigma ** 2)
        )
        g2d /= g2d.max()

        npm, nqm = np.shape(buff)
        spl = _scipy.signal.convolve(buff, g2d)

        return spl[hnbeam : npm + hnbeam, hnbeam : nqm + hnbeam]


class FixedResolutionBufferWhiteNoiseFilter(FixedResolutionBufferFilter):

    """
    This filter adds white noise with the amplitude "bg_lvl" to
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`.
    If "bg_lvl" is not present, 10th percentile of the FRB's value is
    used instead.
    """

    _filter_name = "white_noise"

    def __init__(self, bg_lvl=None):
        self.bg_lvl = bg_lvl

    def apply(self, buff):
        if self.bg_lvl is None:
            amp = np.percentile(buff, 10)
        else:
            amp = self.bg_lvl
        npm, nqm = np.shape(buff)
        return buff + np.random.randn(npm, nqm) * amp
