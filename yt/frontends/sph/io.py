"""
Generic file-handing functions for SPH data




"""
from yt.utilities.io_handler import BaseParticleIOHandler


class IOHandlerSPH(BaseParticleIOHandler):
    """IOHandler implementation specifically for SPH data

    This exists to handle particles with smoothing lengths, which require us
    to read in smoothing lengths along with the the particle coordinates to
    determine particle extents.

    At present this is non-functional.
    """

    pass
