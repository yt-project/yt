from abc import ABC, abstractmethod
from functools import update_wrapper, wraps

import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.visualization.fixed_resolution import FixedResolutionBuffer


def apply_filter(f):
    issue_deprecation_warning(
        "The apply_filter decorator is not used in yt any more and "
        "will be removed in a future version. "
        "Please do not use it.",
        since="4.1",
    )

    @wraps(f)
    def newfunc(self, *args, **kwargs):
        self._filters.append((f.__name__, (args, kwargs)))
        # Invalidate the data of the frb to force its regeneration
        self._data_valid = False
        return self

    return newfunc


class FixedResolutionBufferFilter(ABC):

    """
    This object allows to apply data transformation directly to
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    """

    def __init_subclass__(cls, *args, **kwargs):

        if cls.__init__.__doc__ is None:
            # allow docstring definition at the class level instead of __init__
            cls.__init__.__doc__ = cls.__doc__

        # add a method to FixedResolutionBuffer
        method_name = "apply_" + cls._filter_name

        def closure(self, *args, **kwargs):
            self._filters.append(cls(*args, **kwargs))
            self._data_valid = False
            return self

        update_wrapper(
            wrapper=closure,
            wrapped=cls.__init__,
            assigned=("__annotations__", "__doc__"),
        )

        closure.__name__ = method_name
        setattr(FixedResolutionBuffer, method_name, closure)

    @abstractmethod
    def __init__(self, *args, **kwargs):
        """This method is required in subclasses, but the signature is arbitrary"""
        pass

    @abstractmethod
    def apply(self, buff: np.ndarray) -> np.ndarray:
        pass

    def __call__(self, buff: np.ndarray) -> np.ndarray:
        # alias to apply
        return self.apply(buff)


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
            -((x / sigma) ** 2 + (y / sigma) ** 2) / (2 * sigma**2)
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
