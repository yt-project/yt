from collections import defaultdict

import numpy as np

from yt.frontends.open_pmd.misc import get_component, is_const_component
from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler


class IOHandlerOpenPMDHDF5(BaseIOHandler):
    _field_dtype = "float32"
    _dataset_type = "openPMD"

    def __init__(self, ds, *args, **kwargs):
        self.ds = ds
        self._handle = ds._handle
        self.base_path = ds.base_path
        self.meshes_path = ds.meshes_path
        self.particles_path = ds.particles_path
        self._array_fields = {}
        self._cached_ptype = ""

    def _fill_cache(self, ptype, index=0, offset=None):
        """Fills the particle position cache for the ``ptype``.

        Parameters
        ----------
        ptype : str
            The on-disk name of the particle species
        index : int, optional
        offset : int, optional
        """
        if str((ptype, index, offset)) not in self._cached_ptype:
            self._cached_ptype = str((ptype, index, offset))
            pds = self._handle[self.base_path + self.particles_path + "/" + ptype]
            axes = list(pds["position"].keys())
            if offset is None:
                if is_const_component(pds["position/" + axes[0]]):
                    offset = pds["position/" + axes[0]].attrs["shape"]
                else:
                    offset = pds["position/" + axes[0]].len()
            self.cache = np.empty((3, offset), dtype=np.float64)
            for i in np.arange(3):
                ax = "xyz"[i]
                if ax in axes:
                    np.add(
                        get_component(pds, "position/" + ax, index, offset),
                        get_component(pds, "positionOffset/" + ax, index, offset),
                        self.cache[i],
                    )
                else:
                    # Pad accordingly with zeros to make 1D/2D datasets compatible
                    # These have to be the same shape as the existing axes since that
                    # equals the number of particles
                    self.cache[i] = np.zeros(offset)

    def _read_particle_selection(self, chunks, selector, fields):
        """Read particle fields for particle species masked by a selection.

        Parameters
        ----------
        chunks
            A list of chunks
            A chunk is a list of grids
        selector
            A region (inside your domain) specifying which parts of the field
            you want to read. See [1] and [2]
        fields : array_like
            Tuples (ptype, pfield) representing a field

        Returns
        -------
        dict
            keys are tuples (ptype, pfield) representing a field
            values are (N,) ndarrays with data from that field
        """
        f = self._handle
        bp = self.base_path
        pp = self.particles_path
        ds = f[bp + pp]
        unions = self.ds.particle_unions
        chunks = list(chunks)  # chunks is a generator

        rv = {}
        ind = {}
        particle_count = {}
        ptf = defaultdict(list)  # ParticleTypes&Fields
        rfm = defaultdict(list)  # RequestFieldMapping

        for (ptype, pname) in fields:
            pfield = (ptype, pname)
            # Overestimate the size of all pfields so they include all particles
            # and shrink it later
            particle_count[pfield] = 0
            if ptype in unions:
                for pt in unions[ptype]:
                    particle_count[pfield] += self.ds.particle_type_counts[pt]
                    ptf[pt].append(pname)
                    rfm[pt, pname].append(pfield)
            else:
                particle_count[pfield] = self.ds.particle_type_counts[ptype]
                ptf[ptype].append(pname)
                rfm[pfield].append(pfield)
            rv[pfield] = np.empty((particle_count[pfield],), dtype=np.float64)
            ind[pfield] = 0

        for ptype in ptf:
            for chunk in chunks:
                for grid in chunk.objs:
                    if str(ptype) == "io":
                        species = list(ds.keys())[0]
                    else:
                        species = ptype
                    if species not in grid.ptypes:
                        continue
                    # read particle coords into cache
                    self._fill_cache(species, grid.pindex, grid.poffset)
                    mask = selector.select_points(
                        self.cache[0], self.cache[1], self.cache[2], 0.0
                    )
                    if mask is None:
                        continue
                    pds = ds[species]
                    for field in ptf[ptype]:
                        component = "/".join(field.split("_")[1:])
                        component = component.replace("positionCoarse", "position")
                        component = component.replace("-", "_")
                        data = get_component(pds, component, grid.pindex, grid.poffset)[
                            mask
                        ]
                        for request_field in rfm[(ptype, field)]:
                            rv[request_field][
                                ind[request_field] : ind[request_field] + data.shape[0]
                            ] = data
                            ind[request_field] += data.shape[0]

        for field in fields:
            rv[field] = rv[field][: ind[field]]

        return rv

    def _read_fluid_selection(self, chunks, selector, fields, size):
        """Reads given fields masked by a given selection.

        Parameters
        ----------
        chunks
            A list of chunks
            A chunk is a list of grids
        selector
            A region (inside your domain) specifying which parts of the field
            you want to read. See [1] and [2]
        fields : array_like
            Tuples (fname, ftype) representing a field
        size : int
            Size of the data to read

        Returns
        -------
        dict
            keys are tuples (ftype, fname) representing a field
            values are flat (``size``,) ndarrays with data from that field
        """
        f = self._handle
        bp = self.base_path
        mp = self.meshes_path
        ds = f[bp + mp]
        chunks = list(chunks)

        rv = {}
        ind = {}

        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError

        if size is None:
            size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)
        for field in fields:
            rv[field] = np.empty(size, dtype=np.float64)
            ind[field] = 0

        for (ftype, fname) in fields:
            field = (ftype, fname)
            for chunk in chunks:
                for grid in chunk.objs:
                    mask = grid._get_selector_mask(selector)
                    if mask is None:
                        continue
                    component = fname.replace("_", "/").replace("-", "_")
                    if component.split("/")[0] not in grid.ftypes:
                        data = np.full(grid.ActiveDimensions, 0, dtype=np.float64)
                    else:
                        data = get_component(ds, component, grid.findex, grid.foffset)
                    # The following is a modified AMRGridPatch.select(...)
                    data.shape = (
                        mask.shape
                    )  # Workaround - casts a 2D (x,y) array to 3D (x,y,1)
                    count = grid.count(selector)
                    rv[field][ind[field] : ind[field] + count] = data[mask]
                    ind[field] += count

        for field in fields:
            rv[field] = rv[field][: ind[field]]
            rv[field].flatten()

        return rv
