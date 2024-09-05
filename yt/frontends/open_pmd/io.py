from collections import defaultdict

import numpy as np

from yt.frontends.open_pmd.misc import get_component, is_const_component, coordinate_mapping
from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler


class IOHandlerOpenPMD(BaseIOHandler):
    _field_dtype = "float32"
    _dataset_type = "openPMD"

    def __init__(self, ds, *args, **kwargs):
        self.ds = ds
        self._handle = ds._handle
        self.base_path = ds.base_path
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
            pds = self._handle.particles[ptype]
            axes = list(pds["position"])
            if offset is None:
                offset = pds["position"][axes[0]].shape[0]
            self.cache = np.empty((3, offset), dtype=np.float64)
            for i in np.arange(len(self.ds._axes_labels)): #only in dims of simulation
                #ax = "xyz"[i] # cartesian only
                #this is done to match the transpose and padding done for meshes
                ax = sorted(self.ds._axes_labels)[i]
                #problem, this is where we should correct particle axes order
                #if ax in axes:
                """
                if is_const_component(pds["position"][ax]):
                    # Pad accordingly with 0.5 to make 1D/2D datasets compatible
                    # as missing dimensions are given empty bounds [0,1]
                    # These have to be the same shape as the existing axes since that
                    # equals the number of particles
                    self.cache[i] = np.full(pds["position"][ax].shape[0], 0.5, dtype=np.float64)
                """
                self.cache[i] = np.add(
                    get_component(pds["position"], ax, index, offset),
                    get_component(pds["positionOffset"], ax, index, offset)
                )
            #now the case for dims not in simulation
            for ax in set(['x','y','z']) - set(list(self.ds._axes_labels)):
                # Pad accordingly with 0.5 to make 1D/2D datasets compatible
                # as missing dimensions are given empty bounds [0,1]
                # These have to be the same shape as the existing axes since that
                # equals the number of particles. Only work for Cartesian coordinates
                i+=1
                self.cache[i] = np.full(offset,
                                    0.5, dtype=np.float64)

    #def _read_particle_coords(self, chunks, ptf):
    #    yield from (
    #        (ptype, xyz, 0.0)
    #        for ptype, xyz in self._read_particle_fields(chunks, ptf, None)
    #    )
    
    #def _read_particle_fields(self, chunks, ptf, selector):
    #    

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
        ds = f.particles
        unions = self.ds.particle_unions
        chunks = list(chunks)  # chunks is a generator

        rv = {}
        ind = {}
        particle_count = {}
        ptf = defaultdict(list)  # ParticleTypes&Fields
        rfm = defaultdict(list)  # RequestFieldMapping

        for ptype, pname in fields:
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
            #bit of a hack, but for overlapping 
            rv[pfield] = np.empty((particle_count[pfield],), dtype=np.float64)
            ind[pfield] = 0
            print(ptype, pname, ptf[ptype], rfm[pfield], 'HERENOW', particle_count[pfield])

        for ptype in ptf:
            for chunk in chunks:
                print('len(chunk.objs)', len(chunk.objs))
                for grid in chunk.objs:
                    print(len(chunk.objs), 'new_lengths')
                    if str(ptype) == "io":
                        species = list(ds)[0]
                    else:
                        species = ptype
                    if species not in grid.ptypes:
                        continue
                    # read particle coords into cache
                    #print('pindex,poffset', grid.pindex, grid.poffset, grid.LeftEdge, grid.RightEdge)
                    #grid._read_particles() #fairly inefficient, would be sweet to do better
                    #print('pindex,poffset', grid.pindex, grid.poffset, grid.LeftEdge, grid.RightEdge, grid.RightEdge - grid.LeftEdge)
                    #correctly set our particles indices and offsets according to species. Choice of position is arbitrary here
                    print('lengths', grid.id, len(ds[species]['position'][list(ds[species]['position'])[0]].available_chunks()))
                    grid.pindex  = ds[species]['position'][list(ds[species]['position'])[0]].available_chunks()[grid.id].offset[0]
                    grid.poffset = ds[species]['position'][list(ds[species]['position'])[0]].available_chunks()[grid.id].extent[0]
                    self._fill_cache(species, grid.pindex, grid.poffset) #old

                    #self._fill_cache(species, grid.pindex, grid.poffset) #old
                    print('selector.count', selector.count_points(self.cache[0], self.cache[1], self.cache[2], 0.0))
                    #grid._read_particles()
                    #problem here because don't know what is happening
                    mask = selector.select_points(
                        self.cache[0], #x
                        self.cache[1], #y
                        self.cache[2], #z
                        0.0            #radius
                    )
                    if mask is None:
                        continue
                    grid._read_particles()
                    pds = ds[species]
                    for field in ptf[ptype]:
                        print(field, ptype, "fieldptype")
                        # just chop off 'particle',
                        record = field.split("_")[1:]
                        component = record[-1]
                        if len(record) > 1:  # specifying axes here
                            record = record[-2].replace("positionCoarse", "position")
                            #component = component.replace("-", "_")
                            component = coordinate_mapping(component)
                            #make sure padded dims are in the cell-centered selection region
                            if record == "position" and is_const_component(pds[record][component])\
                            and component not in self.ds._axes_labels:
                                #viewing particles along padded dimensions is fairly useless
                                print("HERE")
                                #data = np.full(pds[record][component].shape[0], 0.5, dtype = np.float64)[grid._pdata[str(species)]['mask']]
                            #else:
                            #print('first5,last5', grid._pdata[str(species)]['mask'][:5], grid._pdata[str(species)]['mask'][-5:])
                            print(pds[record][component].shape[0], 'shape')
                            data = get_component(
                                pds[record],
                                component,
                                grid.pindex,
                                grid.poffset
                            )[mask] #[grid._pdata[str(species)]['mask']]
                            #print('HERE', np.max(data), np.min(data), component)
                        else:  # just particle_mass, charge, weighting,and id
                            print(type(grid.pindex), 'TYPE')
                            data = get_component(
                                pds[component],
                                list(pds[component])[0], #takes care of the SCALAR key for openpmd_api
                                grid.pindex,
                                grid.poffset,
                            )[mask] #[grid._pdata[str(species)]['mask']]
                        #print('mask shapes', np.shape(mask), np.sum(np.diff(grid._pdata[str(species)]['mask'])))
                        for request_field in rfm[ptype, field]:
                            print('shapes', ind[request_field], data.shape, np.shape(data))
                            print(np.shape(rv[request_field]), request_field)
                            rv[request_field][
                                ind[request_field] : ind[request_field] + data.shape[0]
                            ] = data
                            print(np.shape(rv[request_field]))
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
        ds = f.meshes
        chunks = list(chunks)
        rv = {}  # flat fluid array
        ind = {}  # flat indices?
        # this makes me think there will be problems
        #if isinstance(selector, GridSelector):
        #    print()
        #    if not (len(chunks) == len(chunks[1].objs) == 1):
        #        raise RuntimeError
        # print(size, 'presize')
        if size is None:
            size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)
        # print(size)
        for field in fields:
            rv[field] = np.empty(size, dtype=np.float64)
            ind[field] = 0
        for ftype, fname in fields:
            field = (ftype, fname)
            for chunk in chunks:  # one chunk per level?
                for grid in chunk.objs:  # grids per level
                    mask = grid._get_selector_mask(selector)
                    if mask is None:
                        continue
                    if grid.Level > 0:  # we hide grid levels from user
                        component = (
                            fname.split("_")[0]
                            + f"_lvl{str(grid.Level)}_"
                            + fname.split("_")[-1]
                        )
                    else:
                        component = fname
                    component = component.replace("-", "_")
                    if "_".join(fname.split("_")[:-1]) not in grid.ftypes:
                        # we get here due to our last chunk holding just particles
                        data = np.full(grid.ActiveDimensions, 0, dtype=np.float64)
                    else:
                        component_field = "_".join(component.split("_")[:-1])
                        component_axes = component.split("_")[-1]
                        data = get_component(
                            ds[component_field],
                            component_axes,
                            grid.findex.copy(),
                            grid.foffset.copy(),
                        )
                    # The following is a modified AMRGridPatch.select(...)
                    #print(mask.shape, 'mask shape', data.shape, fname, grid.Level, grid.ftypes, grid.ptypes)
                    """
                    THIS WORKS"""
                    #if data.shape[0] != mask.shape[0]:
                    if ds[component_field].data_order == 'C':
                        #this only works for rectangular, non particle grids
                        #print('before swap', np.shape(data))
                        if len(data.shape) == 3:
                            data = np.transpose(data,(2, 1, 0))
                        elif len(data.shape) == 2:
                            data = np.transpose(data.reshape(data.shape[0],data.shape[1],1), (1, 0, 2))
                        elif len(data.shape) == 1:
                            data = data.reshape(data.shape[0],1,1)
                        #data = np.swapaxes(data, 1,0)
                        #data.shape
                        #data = np.transpose(data, _3d)
                        #data = np.swapaxes(data,1,0) #this works when we pad last axes as 1 in 3c
                    #print('datashape,maskshape', data.shape, mask.shape)
                    data.shape = (
                        mask.shape
                    )
                    #print(data.shape, 'postahape', np.shape(data))
                    #OLD VERSION
                    #data.shape = (
                    #   make.shape
                    #)
                    # Workaround - casts a 2D (x,y) array to 3D (x,y,1)
                    #data.reshape(mask.shape)
                    # Now we want (z,x) C orderded daga to be cast to (z,1,x)
                    """
                    Doesn't work"""
                    #data = np.reshape(data, mask.shape)
                    count = grid.count(selector)
                    #print(count,'count')
                    rv[field][ind[field] : ind[field] + count] = data[mask] #multiply by unit SI here?
                    ind[field] += count  # flattened index
        # it would be sweet to flush here
        # how could we do this without looping?
        for field in fields:
            rv[field] = rv[field][: ind[field]]
            rv[field].flatten()
        return rv
 # type: ignore