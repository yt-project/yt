##
## This is a modified version of the LogLocator used in Matplotlib.
## It is subject to the terms of the BSD license, and copyright is held by the
## original authors.
##

import math
import numpy as na

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
            self._subs = na.asarray(subs)+0.0

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
            if numdec>10: subs = na.array([1.0])
            elif numdec>6: subs = na.arange(2.0, b, 2.0)
            else: subs = na.arange(2.0, b)
        else:
            subs = self._subs

        stride = 1
        while numdec/stride+1 > self.numticks:
            stride += 1

        decades = na.arange(math.floor(vmin),
                             math.ceil(vmax)+stride, stride)
        if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
            ticklocs = []
            for decadeStart in b**decades:
                ticklocs.extend( subs*decadeStart )
        else:
            ticklocs = b**decades

        return na.array(ticklocs)

if __name__ == "__main__":
    ll = LogLocator()
    print ll(1e-24, 5e-25)
    print ll(1e-24, 1e-28)
    print ll(1e-24, 1e-35)
