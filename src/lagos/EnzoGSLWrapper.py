"""
Not yet functional, and with SWIG, I'm not sure it ever will be.
"""

from yt.lagos import *
import pygsl

def SolveCoupledChemistryCooling(grid, dt):
    """
    This is a wrapper to the new coupled chemistry cooling module

    @param grid: the EnzoGrid instance to operate upon
    @type grid: L{EnzoGrid<EnzoGrid>}
    @param dt: how long, in code units, to evolve
    @type dt: float
    """
    # First we set up the GSL splining functions
    pass
