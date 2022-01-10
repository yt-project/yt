"""
Generic file-handing functions for SPH data




"""
from typing import DefaultDict, Iterator, List

from yt._typing import SPHParticleCoordinateTuple
from yt.utilities.io_handler import BaseParticleIOHandler


class IOHandlerSPH(BaseParticleIOHandler):
    """IOHandler implementation specifically for SPH data

    This exists to handle particles with smoothing lengths, which require us
    to read in smoothing lengths along with the the particle coordinates to
    determine particle extents.
    """

    def _read_particle_coords(
        self,
        chunks,
        ptf: DefaultDict[str, List[str]],
    ) -> Iterator[SPHParticleCoordinateTuple]:
        # An iterator that yields particle coordinates for each chunk by particle
        # type. Must be implemented by each frontend. Must yield a tuple of
        # (particle type, xyz, smoothing length) by chunk.
        raise NotImplementedError

    def _count_selected_particles(
        self,
        psize: DefaultDict[str, int],
        chunks,
        ptf: DefaultDict[str, List[str]],
        selector,
    ) -> DefaultDict[str, int]:
        # counts the number of particles in a selection by chunk, with smoothing length
        for ptype, (x, y, z), hsml in self._read_particle_coords(chunks, ptf):
            psize[ptype] += selector.count_points(x, y, z, hsml)
        return psize
