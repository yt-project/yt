from yt.frontends.gadget.api import IOHandlerGadgetHDF5
import h5py
import numpy as np

class IOHandlerArepoHDF5(IOHandlerGadgetHDF5):
    _dataset_type = "arepo_hdf5"

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                x = f["/%s/Coordinates" % ptype][si:ei, 0].astype("float64")
                y = f["/%s/Coordinates" % ptype][si:ei, 1].astype("float64")
                z = f["/%s/Coordinates" % ptype][si:ei, 2].astype("float64")
                if ptype == self.ds._sph_ptype:
                    pdtype = f["/%s/Coordinates" % ptype].dtype
                    pshape = f["/%s/Coordinates" % ptype].shape
                    hsml = self._get_smoothing_length(data_file, pdtype, pshape)
                else:
                    hsml = 0.0
                yield ptype, (x, y, z), hsml
            f.close()

    def _get_smoothing_length(self, data_file, position_dtype, position_shape):
        ptype = self.ds._sph_ptype
        ind = int(ptype[-1])
        si, ei = data_file.start, data_file.end
        with h5py.File(data_file.filename, "r") as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][ind].astype("int")
            pcount = np.clip(pcount - si, 0, ei - si)
            ds = f[ptype]["Masses"][si:ei,...]/f[ptype]["Density"][si:ei,...]
            ds *= 3.0/(4.0*np.pi)
            ds **= (1./3.)
            ds *= self.ds.smoothing_factor
            dt = ds.dtype.newbyteorder("N") # Native
            if position_dtype is not None and dt < position_dtype:
                # Sometimes positions are stored in double precision
                # but smoothing lengths are stored in single precision.
                # In these cases upcast smoothing length to double precision
                # to avoid ValueErrors when we pass these arrays to Cython.
                dt = position_dtype
            hsml = np.empty(ds.shape, dtype=dt)
            hsml[:] = ds
            return hsml

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: (x.filename, x.start)):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                g = f["/%s" % ptype]
                if getattr(selector, 'is_all_data', False):
                    mask = slice(None, None, None)
                else:
                    coords = g["Coordinates"][si:ei].astype("float64")
                    if ptype == 'PartType0':
                        hsmls = self._get_smoothing_length(data_file,
                                                           g["Coordinates"].dtype,
                                                           g["Coordinates"].shape)
                    else:
                        hsmls = 0.0
                    mask = selector.select_points(
                            coords[:,0], coords[:,1], coords[:,2], hsmls)
                    del coords
                if mask is None:
                    continue
                for field in field_list:
                    if field in ("Mass", "Masses") and \
                            ptype not in self.var_mass:
                        data = np.empty(mask.sum(), dtype="float64")
                        ind = self._known_ptypes.index(ptype)
                        data[:] = self.ds["Massarr"][ind]
                    elif field in self._element_names:
                        rfield = 'ElementAbundance/' + field
                        data = g[rfield][si:ei][mask, ...]
                    elif field.startswith("Metallicity_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["Metallicity"][si:ei, col][mask]
                    elif field.startswith("Chemistry_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["ChemistryAbundances"][si:ei, col][mask]
                    elif field == "smoothing_length":
                        data = hsmls[mask].astype("float64")
                    else:
                        data = g[field][si:ei][mask, ...]

                    yield (ptype, field), data
            f.close()


    def _identify_fields(self, data_file):
        f = h5py.File(data_file.filename, "r")
        fields = []
        cname = self.ds._particle_coordinates_name  # Coordinates
        mname = self.ds._particle_mass_name  # Mass

        # loop over all keys in OWLS hdf5 file
        #--------------------------------------------------
        for key in f.keys():

            # only want particle data
            #--------------------------------------
            if not key.startswith("PartType"):
                continue

            # particle data group
            #--------------------------------------
            g = f[key]
            if cname not in g:
                continue

            # note str => not unicode!
            ptype = str(key)
            if ptype not in self.var_mass:
                fields.append((ptype, mname))

            # loop over all keys in PartTypeX group
            #----------------------------------------
            for k in g.keys():

                if k == 'ElementAbundance':
                    gp = g[k]
                    for j in gp.keys():
                        kk = j
                        fields.append((ptype, str(kk)))
                elif k == 'Metallicity' and len(g[k].shape) > 1:
                    # Vector of metallicity
                    for i in range(g[k].shape[1]):
                        fields.append((ptype, "Metallicity_%02i" % i))
                elif k == "ChemistryAbundances" and len(g[k].shape) > 1:
                    for i in range(g[k].shape[1]):
                        fields.append((ptype, "Chemistry_%03i" % i))
                else:
                    kk = k
                    if not hasattr(g[kk], "shape"):
                        continue
                    if len(g[kk].shape) > 1:
                        self._vector_fields[kk] = g[kk].shape[1]
                    fields.append((ptype, str(kk)))
        if self._dataset_type == "arepo_hdf5":
            fields.append(("PartType0", "smoothing_length"))
        f.close()
        return fields, {}

