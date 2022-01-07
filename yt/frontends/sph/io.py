"""
Generic file-handing functions for SPH data




"""
from yt.utilities.io_handler import (
    BaseParticleIOHandler,
    ParticleTypeFields,
    ParticleTypeSizes,
)


class IOHandlerSPH(BaseParticleIOHandler):
    """IOHandler implementation specifically for SPH data

    This exists to handle particles with smoothing lengths, which require us
    to read in smoothing lengths along with the the particle coordinates to
    determine particle extents.
    """

    def _count_selected_particles(
        self, psize: ParticleTypeSizes, chunks, ptf: ParticleTypeFields, selector
    ) -> ParticleTypeSizes:
        # counts the number of particles in a selection by chunk, with smoothing length
        for ptype, (x, y, z), hsml in self._read_particle_coords(chunks, ptf):
            psize[ptype] += selector.count_points(x, y, z, hsml)
        return psize
