"""
openPMD-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
import numpy as np


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
        if str((ptype, index, offset)) not in self._cached_ptype:
            # Get a particle species (e.g. /data/3500/particles/e/)
            if ptype in "io":
                spec = self._handle[self.base_path + self.particles_path].keys()[0]
            else:
                spec = ptype
            pds = self._handle[self.base_path + self.particles_path + "/" + spec]
            # Get 1D-Arrays for individual particle positions along axes
            axes = [str(ax) for ax in pds["position"].keys()]
            if self.ds._nonstandard:
                position_offset = "globalCellIdx/"
            else:
                position_offset = "positionOffset/"
            self.cache = np.empty((3, offset), dtype="float64")
            for i in range(3):
                ax = "xyz"[i]
                if ax in axes:
                    np.add(self.get_component(pds, "position/" + ax, index, offset, self.ds._nonstandard),
                           self.get_component(pds, position_offset + ax, index, offset, self.ds._nonstandard),
                           self.cache[i])
                else:
                    # Pad accordingly with zeros to make 1D/2D datasets compatible
                    # These have to be the same shape as the existing axes since that equals the number of particles
                    self.cache[i] = np.zeros(offset)
        self._cached_ptype = str((ptype, index, offset))

    def _read_particle_coords(self, chunks, ptf):
        """
            Reads coordinates for all specified particle-types from file.

            Parameters
            ----------
            chunks:
                A list of chunks
                A chunk is a list of grids

            ptf:
                A dictionary
                - keys are ptypes
                - values are lists of particle fields

            Yields
            ------
            A series of tuples of (ptype_name, (x_coords, y_coords, z_coords)).
            ptype_name is a type pf particle,
            x_coords, y_coords and z_coords are arrays of positions/coordinates of all particles of that type
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
                    for gridptype, idx in g.particle_index:
                        if str(gridptype) == str(spec):
                            index = idx
                    for gridptype, ofs in g.particle_offset:
                        if str(gridptype) == str(spec):
                            offset = ofs
                    mylog.debug("open_pmd - _read_particle_coords: (grid {}) {}, {} [{}:{}]".format(g, ptype, field_list, index, offset))
                    if str((ptype, index, offset)) not in self._cached_ptype:
                        self._fill_cache(ptype, index, offset)
                    yield (ptype, (self.cache[0], self.cache[1], self.cache[2]))

    def _read_particle_fields(self, chunks, ptf, selector):
        """
            Reads particle fields masked by particle type and field out of a list of chunks from file.

            Parameters
            ----------
            chunks:
                A list of chunks
                A chunk is a list of grids

            ptf:
                A dictionary
                - keys are ptype
                - values are lists of (particle?) fields


            selector:
                yt-project.org/docs/dev/quickstart/data_inspection.html?highlight=selector#Examining-Data-in-Regions
                yt-project.org/doc/developing/creating_datatypes.html
                A region (in and/or outside your domain) specifying the field you want to read
                Selector objects have a .select_points(x,y,z) that returns a mask, so you need to do your masking here.

            Yields
            -------
            Tuples of structure ((ptype, field), data)
            where data is a numpy array of data masked to
             - a particle type p
             - in a specified field
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
                    mylog.debug("open_pmd - _read_particle_fields: (grid {}) {}, {} [{}:{}]".format(g, ptype, field_list, index, offset))
                    if str((ptype, index, offset)) not in self._cached_ptype:
                        self._fill_cache(ptype, index, offset)
                    mask = selector.select_points(self.cache[0], self.cache[1], self.cache[2], 0.0)
                    if mask is None:
                        continue
                    pds = ds[self.particles_path + "/" + spec]
                    for field in field_list:
                        nfield = "/".join(field.split("_")[1:]).replace("positionCoarse", "position")
                        data = self.get_component(pds, nfield, index, offset, self.ds._nonstandard)
                        yield ((ptype, field), data[mask])

    def _read_fluid_selection(self, chunks, selector, fields, size):
        """
            Reads selected fields of given size from file.

            From yt doc:
            Receives a collection of data "chunks", a selector describing which "chunks" you are concerned with,
            a list of fields, and the size of the data to read.
            It should create and return a dictionary whose keys are the fields,
            and whose values are numpy arrays containing the data.
            The data should actually be read via the _read_chunk_data() method.

            Parameters
            ----------
            chunks:
                A list of chunks
                A chunk is a list of grids

            selector:
                yt-project.org/docs/dev/quickstart/data_inspection.html?highlight=selector#Examining-Data-in-Regions
                yt-project.org/doc/developing/creating_datatypes.html
                A region (in and/or outside your domain) specifying the field you want to read

            fields:
                A list of (fname, ftype) tuples representing a field

            size:
                Size of the data arrays you want to read

            Returns
            -------
            A dictionary:
            - keys are (ftype, fname) tuples representing a field
            - values are numpy arrays with data form that field
        """
        mylog.info("_read_fluid_selection {} {} {} {}".format(chunks, selector, fields, size))
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
                    nfield = fname.replace("_", "/").replace("-","_")
                    index = g.mesh_index
                    offset = g.mesh_offset
                    data = np.array(self.get_component(ds, nfield, index, offset, self.ds._nonstandard))
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
        """Grab a component from a group

        Parameters
        ----------
        group: hdf5 group to get whole/masked component from
        component_name: string of a component (relative path) inside the group
        index: (optional) start index of data to return
        offset: (optional) size of the data to return
        nonstandard: (optional) bool to signal whether to assume (non)standard PMD

        Returns
        -------
        n-dimensional numpy array filled with values of component
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
        if "value" in record_component.attrs.keys():
            if offset is not None:
                shape = offset
            else:
                shape = record_component.attrs["shape"] - index
            # component is constant, craft an array by hand
            mylog.debug("open_pmd - misc - get_component (const): {}/{}({})".format(group.name, component_name, shape))
            return np.full(shape, record_component.attrs["value"] * unitSI)
        else:
            if offset is not None:
                offset += index
            # component is a dataset, return it (possibly masked)
            mylog.debug(
                "open_pmd - misc - get_component: {}/{}[{}:{}]".format(group.name, component_name, index, offset))
            return np.multiply(record_component[index:offset], unitSI)