"""
openPMD-specific IO functions



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np

from yt.frontends.open_pmd.misc import is_const_component
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


class IOHandlerOpenPMD(BaseIOHandler):
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
        index : int, optional
        offset : int, optional
        """
        if str((ptype, index, offset)) not in self._cached_ptype:
            if ptype in "io":
                spec = self._handle[self.base_path + self.particles_path].keys()[0]
            else:
                spec = ptype
            pds = self._handle[self.base_path + self.particles_path + "/" + spec]
            axes = [str(ax) for ax in pds["position"].keys()]
            if self.ds._nonstandard:
                position_offset = "globalCellIdx/"
            else:
                position_offset = "positionOffset/"
            if offset is None:
                if is_const_component(pds["position/" + axes[0]]):
                    size = pds["position/" + axes[0]].attrs["shape"]
                else:
                    size = pds["position/" + axes[0]].len()
            else:
                size = offset
            self.cache = np.empty((3, size), dtype="float64")
            for i in range(3):
                ax = "xyz"[i]
                if ax in axes:
                    np.add(self.get_component(pds, "position/" + ax, index, offset, self.ds._nonstandard),
                           self.get_component(pds, position_offset + ax, index, offset, self.ds._nonstandard),
                           self.cache[i])
                else:
                    # Pad accordingly with zeros to make 1D/2D datasets compatible
                    # These have to be the same shape as the existing axes since that equals the number of particles
                    self.cache[i] = np.zeros(size)
        self._cached_ptype = str((ptype, index, offset))

    def _read_particle_coords(self, chunks, ptf):
        """Reads coordinates for given particle-types in given chunks from file.

        Parameters
        ----------
        chunks
            A list of chunks
            A chunk is a list of grids
        ptf : dict
            keys are ptypes
            values are lists of particle fields

        Yields
        ------
        tuple : (str, ((N,) ndarray, (N,) ndarray, (N,) ndarray))
            Tuple of ptype and tuple of coordinates (x, y, z) corresponding to coordinates of particles of that ptype
            All coordinate arrays have the same length
        """

        chunks = list(chunks)
        f = self._handle
        ds = f[self.base_path]
        for chunk in chunks:
            for g in chunk.objs:
                if g.filename is None:
                    continue
                for ptype, field_list in sorted(ptf.items()):
                    if ptype in "io":
                        spec = ds[self.particles_path].keys()[0]
                    else:
                        spec = ptype
                    index = dict(g.particle_index).get(spec)
                    offset = dict(g.particle_offset).get(spec)
                    mylog.debug(
                        "open_pmd - _read_particle_coords: (grid {}) {}, {} [{}:{}]".format(g, ptype, field_list, index,
                                                                                            offset))
                    if str((ptype, index, offset)) not in self._cached_ptype:
                        self._fill_cache(ptype, index, offset)
                    yield (ptype, (self.cache[0], self.cache[1], self.cache[2]))

    def _read_particle_fields(self, chunks, ptf, selector):
        """Reads given fields for given particle types masked by a given selection.

        Parameters
        ----------
        chunks
            A list of chunks
            A chunk is a list of grids
        ptf : dict
            keys are ptype
            values are lists of particle fields
        selector
            A region (inside your domain) specifying which parts of the field you want to read
            See [1] and [2]

        References
        ----------
        .. [1] yt-project.org/docs/dev/quickstart/data_inspection.html?highlight=selector#Examining-Data-in-Regions
        .. [2] yt-project.org/doc/developing/creating_datatypes.html

        Yields
        ------
        tuple : ((str, str), (N,) ndarray)
            Tuple of tuple (ptype, fieldname) and masked field-data
        """
        chunks = list(chunks)
        f = self._handle
        ds = f[self.base_path]
        for chunk in chunks:
            for g in chunk.objs:
                if g.filename is None:
                    continue
                for ptype, field_list in sorted(ptf.items()):
                    # Get a particle species (e.g. /data/3500/particles/e/)
                    if ptype in "io":
                        spec = ds[self.particles_path].keys()[0]
                    else:
                        spec = ptype
                    index = dict(g.particle_index).get(spec)
                    offset = dict(g.particle_offset).get(spec)
                    mylog.debug(
                        "open_pmd - _read_particle_fields: (grid {}) {}, {} [{}:{}]".format(g, ptype, field_list, index,
                                                                                            offset))
                    if str((ptype, index, offset)) not in self._cached_ptype:
                        self._fill_cache(ptype, index, offset)
                    mask = selector.select_points(self.cache[0], self.cache[1], self.cache[2], 0.0)
                    if mask is None:
                        continue
                    pds = ds[self.particles_path + "/" + spec]
                    for field in field_list:
                        component = "/".join(field.split("_")[1:]).replace("positionCoarse",
                                                                           "position").replace("-", "_")
                        data = self.get_component(pds, component, index, offset, self.ds._nonstandard)
                        yield ((ptype, field), data[mask])

    def _read_fluid_selection(self, chunks, selector, fields, size):
        """Reads given fields for given meshes masked by a given selection.

        Parameters
        ----------
        chunks
            A list of chunks
            A chunk is a list of grids
        selector
            A region (inside your domain) specifying which parts of the field you want to read
            See [1] and [2]
        fields : array_like
            tuples (fname, ftype) representing a field
        size : int
            Size of the data to read

        Returns
        -------
        dict
            keys are tuples (ftype, fname) representing a field
            values are flat (``size``,) ndarrays with data from that field
        """
        mylog.debug("open_pmd - _read_fluid_selection {} {} {} {}".format(chunks, selector, fields, size))
        f = self._handle
        bp = self.base_path
        mp = self.meshes_path
        rv = {}
        chunks = list(chunks)

        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError

        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        ind = {}
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")
            ind[field] = 0

        for ftype, fname in fields:
            field = (ftype, fname)
            for chunk in chunks:
                for g in chunk.objs:
                    ds = f[bp + mp]
                    component = fname.replace("_", "/").replace("-", "_")
                    index = g.mesh_index
                    offset = g.mesh_offset
                    data = np.array(self.get_component(ds, component, index, offset, self.ds._nonstandard))
                    # The following is a modified AMRGridPatch.select(...)
                    mask = g._get_selector_mask(selector)
                    data.shape = mask.shape  # Workaround - casts a 2D (x,y) array to 3D (x,y,1)
                    count = g.count(selector)
                    rv[field][ind[field]:ind[field] + count] = data[mask]
                    ind[field] += count
        for i in rv:
            rv[i].flatten()
        return rv

    def get_component(self, group, component_name, index=0, offset=None, nonstandard=False):
        """Grab a dataset component from a group as a whole or sliced.

        Parameters
        ----------
        group : h5py.Group
        component_name : str
            relative path of the component in the group
        index : int, optional
            first entry along the first axis to read
        offset : int, optional
            number of entries to read
            if not supplied, every entry after index is returned
        nonstandard : bool, optional

        Notes
        -----
        This scales every entry of the component with the respective "unitSI".

        Returns
        -------
        (N,) ndarray

        """
        record_component = group[component_name]
        if nonstandard:
            try:
                unitSI = record_component.attrs["sim_unit"]
                if "globalCellIdx" in component_name:
                    it = group["/data"].keys()[0]
                    length = group["/data/" + it].attrs['unit_length']
                    unitSI *= length
            except:
                unitSI = 1.0
        else:
            unitSI = record_component.attrs["unitSI"]
        # check whether component is constant
        if is_const_component(record_component):
            if offset is not None:
                shape = offset
            else:
                shape = record_component.attrs["shape"] - index
            # component is constant, craft an array by hand
            mylog.debug("open_pmd - get_component (const): {}/{}({})".format(group.name, component_name, shape))
            return np.full(shape, record_component.attrs["value"] * unitSI)
        else:
            if offset is not None:
                offset += index
            # component is a dataset, return it (possibly masked)
            mylog.debug(
                "open_pmd - get_component: {}/{}[{}:{}]".format(group.name, component_name, index, offset))
            return np.multiply(record_component[index:offset], unitSI)
