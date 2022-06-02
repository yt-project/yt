"""
Generic file-handing functions for SPH data




"""
from yt.utilities.io_handler import BaseParticleIOHandler


class IOHandlerSPH(BaseParticleIOHandler):
    """IOHandler implementation specifically for SPH data

    This exists to handle particles with smoothing lengths, which require us
    to read in smoothing lengths along with the the particle coordinates to
    determine particle extents.
    """

    def _get_smoothing_length(
        self, data_file, position_dtype, position_shape, handle=None
    ):
        raise NotImplementedError("SPH frontends must implement _get_smoothing_length")
