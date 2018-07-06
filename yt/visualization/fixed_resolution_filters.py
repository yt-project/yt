"""
Fixed resolution buffer filters.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.extern.six import add_metaclass
from functools import wraps

filter_registry = {}


def apply_filter(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        frb = args[0]
        frb._filters.append((f.__name__, (args, kwargs)))
        # Invalidate the data of the frb to force its regeneration
        frb._data_valid = False
        return args[0]
    return newfunc


class RegisteredFilter(type):

    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        filter_registry[name] = cls


@add_metaclass(RegisteredFilter)
class FixedResolutionBufferFilter(object):

    """
    This object allows to apply data transformation directly to
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`
    """

    def __init__(self, *args, **kwargs):
        pass

    def apply(self, buff):
        pass


class FixedResolutionBufferGaussBeamFilter(FixedResolutionBufferFilter):

    """
    This filter convolves
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer` with
    2d gaussian that has standard deviation 'sigma'. Extra arguments
    are directly passed to `scipy.ndimage.gaussian_filter`.
    """
    _filter_name = 'gauss_beam'

    def __init__(self, sigma=2.0, **kwa):
        self.sigma = sigma
        self.extra_args = kwa

    def apply(self, buff):
        from yt.utilities.on_demand_imports import _scipy
        spl = _scipy.ndimage.gaussian_filter(buff, self.sigma, **self.extra_args)
        return spl


class FixedResolutionBufferWhiteNoiseFilter(FixedResolutionBufferFilter):

    """
    This filter adds white noise with the amplitude "bg_lvl" to
    :class:`yt.visualization.fixed_resolution.FixedResolutionBuffer`.
    If "bg_lvl" is not present, 10th percentile of the FRB's value is
    used instead.
    """
    _filter_name = 'white_noise'

    def __init__(self, bg_lvl=None):
        self.bg_lvl = bg_lvl

    def apply(self, buff):
        if self.bg_lvl is None:
            amp = np.percentile(buff, 10)
        else:
            amp = self.bg_lvl
        npm, nqm = np.shape(buff)
        return buff + np.random.randn(npm, nqm) * amp
