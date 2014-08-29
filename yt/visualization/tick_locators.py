from __future__ import print_function
##
## This is a modified version of the LogLocator used in Matplotlib.
## It is subject to the terms of the BSD license, and copyright is held by the
## original authors.
##

import math
import numpy as np

def is_decade(x,base=10):
    if x == 0.0:
        return True
    lx = math.log(x)/math.log(base)
    return lx==int(lx)

class LogLocator(object):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, subs=[1.0], numdecs=4):
        """
        place ticks on the location= base**i*subs[j]
        """
        self.base(base)
        self.subs(subs)
        self.numticks = 15
        self.numdecs = numdecs

    def base(self,base):
        """
        set the base of the log scaling (major tick every base**i, i interger)
        """
        self._base=base+0.0

    def subs(self,subs):
        """
        set the minor ticks the log scaling every base**i*subs[j]
        """
        if subs is None:
            self._subs = None  # autosub
        else:
            self._subs = np.asarray(subs)+0.0

    def _set_numticks(self):
        self.numticks = 15  # todo; be smart here; this is just for dev

    def __call__(self, vmin, vmax):
        'Return the locations of the ticks'
        b=self._base

        if vmin <= 0.0:
            raise ValueError(
                "Data has no positive values, and therefore can not be log-scaled.")

        vmin = math.log(vmin)/math.log(b)
        vmax = math.log(vmax)/math.log(b)

        if vmax<vmin:
            vmin, vmax = vmax, vmin

        numdec = math.floor(vmax)-math.ceil(vmin)

        if self._subs is None: # autosub
            if numdec>10: subs = np.array([1.0])
            elif numdec>6: subs = np.arange(2.0, b, 2.0)
            else: subs = np.arange(2.0, b)
        else:
            subs = self._subs

        stride = 1
        while numdec/stride+1 > self.numticks:
            stride += 1

        decades = np.arange(math.floor(vmin),
                             math.ceil(vmax)+stride, stride)
        if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
            ticklocs = []
            for decadeStart in b**decades:
                ticklocs.extend( subs*decadeStart )
        else:
            ticklocs = b**decades

        return np.array(ticklocs)


class LinearLocator(object):
    """
    Determine the tick locations

    The first time this function is called it will try to set the
    number of ticks to make a nice tick partitioning.  Thereafter the
    number of ticks will be fixed so that interactive navigation will
    be nice
    """


    def __init__(self, numticks = None, presets=None):
        """
        Use presets to set locs based on lom.  A dict mapping vmin, vmax->locs
        """
        self.numticks = numticks
        if presets is None:
            self.presets = {}
        else:
            self.presets = presets

    def __call__(self, vmin, vmax):
        'Return the locations of the ticks'

        # vmin, vmax = self.axis.get_view_interval()
        # vmin, vmax = mtransforms.nonsingular(vmin, vmax, expander = 0.05)
        if vmax<vmin:
            vmin, vmax = vmax, vmin

        if (vmin, vmax) in self.presets:
            return self.presets[(vmin, vmax)]

        if self.numticks is None:
            self._set_numticks()



        if self.numticks==0: return []
        ticklocs = np.linspace(vmin, vmax, self.numticks)

        #return self.raise_if_exceeds(ticklocs)
        return ticklocs


    def _set_numticks(self):
        self.numticks = 11  # todo; be smart here; this is just for dev

    # def view_limits(self, vmin, vmax):
    #     'Try to choose the view limits intelligently'

    #     if vmax<vmin:
    #         vmin, vmax = vmax, vmin

    #     if vmin==vmax:
    #         vmin-=1
    #         vmax+=1

    #     exponent, remainder = divmod(math.log10(vmax - vmin), 1)

    #     if remainder < 0.5:
    #         exponent -= 1
    #     scale = 10**(-exponent)
    #     vmin = math.floor(scale*vmin)/scale
    #     vmax = math.ceil(scale*vmax)/scale

    #     return mtransforms.nonsingular(vmin, vmax)


if __name__ == "__main__":
    ll = LogLocator()
    print(ll(1e-24, 5e-25))
    print(ll(1e-24, 1e-28))
    print(ll(1e-24, 1e-35))
    lll = LinearLocator()
    print(lll(-1e-24, 1e-24))
    print(lll(-2.3, 1.3))
    print(lll(10,23.))
