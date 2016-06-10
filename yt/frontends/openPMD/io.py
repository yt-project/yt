"""
openPMD-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
# Copyright (c) 2015, Daniel Grassinger (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.extern.six import u, b, iteritems
from yt.geometry.selection_routines import mask_fill, AlwaysSelector
from yt.utilities.logger import ytLogger as mylog
import h5py
import numpy as np
from collections import defaultdict
from .data_structures import openPMDBasePath
from .misc import *


class IOHandlerOpenPMD(BaseIOHandler, openPMDBasePath):
    # TODO Data should be loaded chunk-wise to support parallelism. Fields can be chunked arbitrarily in space, particles should be read via particlePatches.
    # TODO Maybe reuse hdf5 file handles in loops?
    # TODO Change strategies when dealing with different iterationEncoding?
    # TODO since the standard allows \w+ regexes, replace _'s in oPMD for yt

    # TODO Unify numpy data types?
    # THERE WAS EVEN A GLOBAL VARIABLE FOR THIS
    _field_dtype = "float32"
    _dataset_type = "openPMD"

    def __init__(self, ds, *args, **kwargs):

        self.ds = ds
        self._handle = ds._handle
        self._setBasePath(self._handle, self.ds._filepath)
        self.meshesPath = self._handle["/"].attrs["meshesPath"]
        self.particlesPath = self._handle["/"].attrs["particlesPath"]
        self._cached_ptype = ""
        self._cache = {}
        self._component_getter = ComponentGetter()

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
        self._array_fields = {}
        for chunk in chunks:
            for g in chunk.objs:
                if g.filename is None:
                    continue
                f = self._handle

                dds = f[self.basePath]
                for ptype, field_list in sorted(ptf.items()):
                    mylog.debug("openPMD - _read_particle_coords: {}, {}".format(ptype, field_list))
                    if ptype in self._cached_ptype:
                        # use the cached x,y,z
                        yield (ptype, (self._cache['x'], self._cache['y'], self._cache['z']))
                    else:
                        # Get a particle species (e.g. /data/3500/particles/e/)
                        if "io" in ptype:
                            spec = dds[self.particlesPath].keys()[0]
                        else:
                            spec = ptype
                        pds = dds[self.particlesPath + "/" + spec]

                        # Get single 1D-arrays of coordinates from all particular particles
                        axes = [str(ax) for ax in pds["position"].keys()]
                        for i in 'xyz':
                            if i not in axes:
                                # TODO Get 2D particle data working
                                mylog.error("openPMD - Can not handle 2D particle data for the moment!")
                                continue
                        pos = {}
                        off = {}
                        for ax in axes:
                            pos[ax] = self._component_getter.get_component(pds, "position/" + ax)
                            off[ax] = self._component_getter.get_component(pds, "positionOffset/" + ax)
                            self._cache[ax] = pos[ax] + off[ax]
                        self._cached_ptype = ptype

                        for field in field_list:
                            # Trim field path (yt scheme) to match openPMD scheme
                            nfield = "/".join(field.split("_")[1:])

                            # Save the size of the particle fields
                            # TODO This CAN be used for speedup, not sure how
                            if is_const_component(pds[nfield]):
                                shape = np.asarray(pds[nfield].attrs["shape"])
                            else:
                                shape = np.asarray(pds[nfield].shape)
                            if shape.ndim > 1:
                                self._array_fields[field] = shape
                        yield (ptype, (self._cache['x'], self._cache['y'], self._cache['z']))

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
        for chunk in chunks:
            for g in chunk.objs:
                if g.filename is None:
                    continue
                f = self._handle

                ds = f[self.basePath]
                for ptype, field_list in sorted(ptf.items()):
                    mylog.debug("openPMD - _read_particle_fields: {}, {}".format(ptype, field_list))

                    # Get a particle species (e.g. /data/3500/particles/e/)
                    if "io" in ptype:
                        spec = ds[self.particlesPath].keys()[0]
                    else:
                        spec = ptype
                    pds = ds[self.particlesPath + "/" + spec]

                    if ptype not in self._cached_ptype:
                        # Get single 1D-arrays of coordinates from all particular particles
                        axes = [str(ax) for ax in pds["position"].keys()]
                        for i in 'xyz':
                            if i not in axes:
                                # TODO Get 2D particle data working
                                mylog.error("openPMD - Can not handle 2D particle data for the moment!")
                                continue
                        pos = {}
                        off = {}
                        for ax in axes:
                            pos[ax] = self._component_getter.get_component(pds, "position/" + ax)
                            off[ax] = self._component_getter.get_component(pds, "positionOffset/" + ax)
                            self._cache[ax] = pos[ax] + off[ax]
                        self._cached_ptype = ptype

                    mask = selector.select_points(self._cache['x'], self._cache['y'], self._cache['z'], 0.0)
                    if mask is None:
                        continue
                    for field in field_list:
                        nfield = "/".join(field.split("_")[1:])
                        data = self._component_getter.get_component(pds, nfield)
                        yield ((ptype, field), data[mask])


    # def _read_particle_selection(self, chunks, selector, fields):
    #     """
    #         Reads masked selection of particles from specified fields.
    #
    #         Receives a collection of data "chunks", a selector describing which "chunks" you are concerned with and
    #         a list of fields.
    #         It should create and return a dictionary whose keys are the fields,
    #         and whose values are numpy arrays containing the data.
    #         The data should actually be read via the _read_chunk_data() method.
    #
    #         Parameters
    #         ----------
    #         chunks:
    #             A list of chunks
    #             A chunk is a list of grids
    #
    #         selector:
    #             yt-project.org/docs/dev/quickstart/data_inspection.html?highlight=selector#Examining-Data-in-Regions
    #             yt-project.org/doc/developing/creating_datatypes.html
    #             A region (in and/or outside your domain) specifying the field you want to read
    #
    #         fields:
    #             A list of (fname, ftype) tuples representing a field
    #
    #         size:
    #             Size of the data arrays you want to read
    #
    #         Returns
    #         -------
    #         A dictionary:
    #         - keys are (ftype, fname) tuples representing a field
    #         - values are numpy arrays with data form that field
    #     """
    #     mylog.debug("Read particle selection")
    #     rv = {}
    #     chunks = list(chunks)
    #
    #     if selector.__class__.__name__ == "GridSelector":
    #         if not (len(chunks) == len(chunks[0].objs) == 1):
    #             raise RuntimeError
    #
    #     rv = {f: np.array([]) for f in fields}
    #     for field in fields:
    #         ftype, fname = field
    #         for chunk in chunks:
    #             for g in chunk.objs:
    #                 mylog.debug("g.select({}, {}, rv[{}], {})".format(selector, data, field,ind))
    #                 rv.update(self._read_chunk_data(chunk, fields))
    #     return rv
    #     # for chunk in chunks:
    #     #     for i in self._read_chunk_data(chunk,fields):
    #     #         mylog.info("self._read_chunk_data({},{}):{}".format(chunk,fields,i))
    #
    #         # for grid in chunk.objs:
    #         #     for ftype, fname in fields:
    #         #         data = self._read_data((ftype, fname))
    #         #         rv[ftype, fname] = np.concatenate((data, rv[ftype, fname]))
    #
    #     """
    #     From Exo 2
    #     for field in fields:
    #             for g in chunk.objs:
    #                             ind += g.select(selector, data, rv[field], ind)
    #     """

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
        mylog.info("_read_fluid_selection")
        rv = {}
        chunks = list(chunks)
        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            # TODO We can probably skip this since the following loops do this, anyway
            # g = chunks[0].objs[0]
            # rv.update(self._read_chunk_data(chunks[0],fields))
            #f = h5py.File(u(g.filename), 'r')
            #for ftype, fname in fields:
            #    rv[ftype, fname] = self._read_data(g, fname)
            #f.close()
            # return rv



        if size is None:
            size = sum((g.count(selector) for chunk in chunks
                        for g in chunk.objs))
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")

        for chunk in chunks:
            rv.update(self._read_chunk_data(chunk, fields))
            # for g in chunk.objs:
                # # WHY do we create a filehandle but not use it?
                # f = h5py.File(u(g.filename), 'r')
                # ###self._handle = f
                # for ftype, fname in fields:
                #     # WHY call with g (a grid) as self?
                #     rv[ftype, fname] = self._read_data(g, fname)
                # f.close()
        return rv

    def _read_data(self, field):
        """
            Reads data from file belonging to a (mesh or particle) field

            Parameters
            ----------
            fieldname:
                Field to get the data for

            Returns
            -------
            A flat numpy array
        """
        ftype, fname = field
        mylog.debug("openPMD - _read_data - reading field: {}.{}".format(ftype, fname))

        if ftype in self.ds.particle_types:
            if "io" in ftype:
                # There is only one particle species, use that one ("io" does not suffice)
                spec = self._handle[self.basePath + self.particlesPath].keys()[0]
            else:
                # Use the provided ftype that corresponds to the species name in the hdf5 file
                spec = str(ftype)
            ds = self._handle[self.basePath + self.particlesPath].get(spec)
            nfield = "/".join(fname.split("_")[1:]).replace("-","_")
        else:
            ds = self._handle[self.basePath + self.meshesPath]
            nfield = fname.replace("_", "/").replace("-","_")
        data = self._component_getter.get_component(ds, nfield)
        return np.array(data).flatten()

    def _read_chunk_data(self, chunk, fields):
        """
            Reads fields specified in the chunk from disk.

            From yt doc:
            Receives a "chunk" of data along with a list of fields we want to read. It loops over all the grid objects
            within the "chunk" of data and reads from disk the specific fields,
            returning a dictionary whose keys are the fields and whose values are numpy arrays of the data.

            Parameters
            ----------
            chunk:
                A single chunk is a list of grids

            fields:
                A list of fields to be read

            Returns
            -------
            A dictionary:
            - keys are (ftype, fname) field tuples
            - values are flat numpy arrays with data form that field
        """
        mylog.debug("_read_chunk_data")
        rv = {}
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append((ftype, fname))
            else:
                fluid_fields.append((ftype, fname))
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update(self._read_particle_selection(
                [chunk], selector, particle_fields))
        if len(fluid_fields) == 0: return rv
        grids_by_file = defaultdict(list)
        # Read fluid fields
        for g in chunk.objs:
            if g.filename is None:
                continue
            grids_by_file[g.filename].append(g)
        for filename in grids_by_file:
            grids = grids_by_file[filename]
            grids.sort()
            f = h5py.File(filename, 'r')
            for g in grids:
                for ftype, fname in fluid_fields:
                    # TODO update basePath for different files
                    data = self._read_data((ftype, fname))
                    # data = f[self.basePath + self.meshPath + fname.replace("_", "/").replace("-","_")]
                    rv[(ftype,fname)] = np.array(data).flatten()
            f.close()
        return rv
