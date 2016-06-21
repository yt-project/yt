"""
OWLS data-file handling function




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton

from .definitions import \
    ghdf5_ptypes

CHUNKSIZE = 10000000

def _get_h5_handle(fn):
    try:
        f = h5py.File(fn, "r")
    except IOError:
        print("ERROR OPENING %s" % (fn))
        if os.path.exists(fn):
            print("FILENAME EXISTS")
        else:
            print("FILENAME DOES NOT EXIST")
        raise
    return f

class IOHandlerOWLS(BaseIOHandler):
    _dataset_type = "OWLS"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")
    _known_ptypes = ghdf5_ptypes
    _var_mass = None
    _element_names = ('Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen',
                       'Neon', 'Magnesium', 'Silicon', 'Iron' )


    @property
    def var_mass(self):
        if self._var_mass is None:
            vm = []
            for i, v in enumerate(self.ds["Massarr"]):
                if v == 0:
                    vm.append(self._known_ptypes[i])
            self._var_mass = tuple(vm)
        return self._var_mass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: x.filename):
            f = _get_h5_handle(data_file.filename)
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                x = f["/%s/Coordinates" % ptype][:,0].astype("float64")
                y = f["/%s/Coordinates" % ptype][:,1].astype("float64")
                z = f["/%s/Coordinates" % ptype][:,2].astype("float64")
                yield ptype, (x, y, z)
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: x.filename):
            f = _get_h5_handle(data_file.filename)
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                g = f["/%s" % ptype]
                coords = g["Coordinates"][:].astype("float64")
                mask = selector.select_points(
                            coords[:,0], coords[:,1], coords[:,2], 0.0)
                del coords
                if mask is None: continue
                for field in field_list:

                    if field in ("Mass", "Masses") and \
                        ptype not in self.var_mass:
                        data = np.empty(mask.sum(), dtype="float64")
                        ind = self._known_ptypes.index(ptype)
                        data[:] = self.ds["Massarr"][ind]

                    elif field in self._element_names:
                        rfield = 'ElementAbundance/' + field
                        data = g[rfield][:][mask,...]
                    elif field.startswith("Metallicity_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["Metallicity"][:,col][mask]
                    elif field.startswith("Chemistry_"):
                        col = int(field.rsplit("_", 1)[-1])
                        data = g["ChemistryAbundances"][:,col][mask]
                    else:
                        data = g[field][:][mask,...]

                    yield (ptype, field), data
            f.close()

    def _initialize_index(self, data_file, regions):
        index_ptype = self.index_ptype
        f = _get_h5_handle(data_file.filename)
        if index_ptype == "all":
            pcount = f["/Header"].attrs["NumPart_ThisFile"][:].sum()
            keys = f.keys()
        else:
            pt = int(index_ptype[-1])
            pcount = f["/Header"].attrs["NumPart_ThisFile"][pt]
            keys = [index_ptype]
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        for key in keys:
            if not key.startswith("PartType"): continue
            if "Coordinates" not in f[key]: continue
            ds = f[key]["Coordinates"]
            dt = ds.dtype.newbyteorder("N") # Native
            pos = np.empty(ds.shape, dtype=dt)
            pos[:] = ds
            regions.add_data_file(pos, data_file.file_id,
                                  data_file.ds.filter_bbox)
            morton[ind:ind+pos.shape[0]] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge,
                data_file.ds.filter_bbox)
            ind += pos.shape[0]
        f.close()
        return morton

    def _count_particles(self, data_file):
        f = _get_h5_handle(data_file.filename)
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:]
        f.close()
        npart = dict(("PartType%s" % (i), v) for i, v in enumerate(pcount))
        return npart


    def _identify_fields(self, data_file):
        f = _get_h5_handle(data_file.filename)
        fields = []
        cname = self.ds._particle_coordinates_name  # Coordinates
        mname = self.ds._particle_mass_name  # Mass

        # loop over all keys in OWLS hdf5 file
        #--------------------------------------------------
        for key in f.keys():

            # only want particle data
            #--------------------------------------
            if not key.startswith("PartType"): continue

            # particle data group
            #--------------------------------------
            g = f[key]
            if cname not in g: continue

            # note str => not unicode!

            #ptype = int(key[8:])
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
                elif k == "ChemistryAbundances" and len(g[k].shape)>1:
                    for i in range(g[k].shape[1]):
                        fields.append((ptype, "Chemistry_%03i" % i))
                else:
                    kk = k
                    if not hasattr(g[kk], "shape"): continue
                    fields.append((ptype, str(kk)))


        f.close()
        return fields, {}
