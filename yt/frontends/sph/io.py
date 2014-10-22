"""
Gadget-specific data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import h5py
import numpy as np
from .definitions import gadget_ptypes, ghdf5_ptypes
from yt.funcs import *
from yt.utilities.exceptions import *

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.fortran_utils import read_record
from yt.utilities.lib.geometry_utils import compute_morton

from yt.geometry.oct_container import _ORDER_MAX

try:
    import requests
except ImportError:
    requests = None

CHUNKSIZE = 10000000

def _get_h5_handle(fn):
    try:
        f = h5py.File(fn, "r")
    except IOError as e:
        print "ERROR OPENING %s" % (fn)
        if os.path.exists(fn):
            print "FILENAME EXISTS"
        else:
            print "FILENAME DOES NOT EXIST"
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
        for data_file in sorted(data_files):
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
        for data_file in sorted(data_files):
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
        f = _get_h5_handle(data_file.filename)
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].sum()
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        for key in f.keys():
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

class IOHandlerEagleNetwork(IOHandlerOWLS):
    _dataset_type = "eagle_network"

class IOHandlerGadgetHDF5(IOHandlerOWLS):
    _dataset_type = "gadget_hdf5"

ZeroMass = object()

class IOHandlerGadgetBinary(BaseIOHandler):
    _dataset_type = "gadget_binary"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")

    # Particle types (Table 3 in GADGET-2 user guide)
    #
    # Blocks in the file:
    #   HEAD
    #   POS
    #   VEL
    #   ID
    #   MASS    (variable mass only)
    #   U       (gas only)
    #   RHO     (gas only)
    #   HSML    (gas only)
    #   POT     (only if enabled in makefile)
    #   ACCE    (only if enabled in makefile)
    #   ENDT    (only if enabled in makefile)
    #   TSTP    (only if enabled in makefile)

    _var_mass = None

    def __init__(self, ds, *args, **kwargs):
        self._fields = ds._field_spec
        self._ptypes = ds._ptype_spec
        super(IOHandlerGadgetBinary, self).__init__(ds, *args, **kwargs)

    @property
    def var_mass(self):
        if self._var_mass is None:
            vm = []
            for i, v in enumerate(self.ds["Massarr"]):
                if v == 0:
                    vm.append(self._ptypes[i])
            self._var_mass = tuple(vm)
        return self._var_mass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype in ptf:
                # This is where we could implement sub-chunking
                f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                pos = self._read_field_from_file(f,
                            tp[ptype], "Coordinates")
                yield ptype, (pos[:,0], pos[:,1], pos[:,2])
            f.close()

    def _read_particle_fields(self, chunks, ptf, selector):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items()):
                f.seek(poff[ptype, "Coordinates"], os.SEEK_SET)
                pos = self._read_field_from_file(f,
                            tp[ptype], "Coordinates")
                mask = selector.select_points(
                    pos[:,0], pos[:,1], pos[:,2], 0.0)
                del pos
                if mask is None: continue
                for field in field_list:
                    if field == "Mass" and ptype not in self.var_mass:
                        data = np.empty(mask.sum(), dtype="float64")
                        m = self.ds.parameters["Massarr"][
                            self._ptypes.index(ptype)]
                        data[:] = m
                        yield (ptype, field), data
                        continue
                    f.seek(poff[ptype, field], os.SEEK_SET)
                    data = self._read_field_from_file(f, tp[ptype], field)
                    data = data[mask,...]
                    yield (ptype, field), data
            f.close()

    def _read_field_from_file(self, f, count, name):
        if count == 0: return
        if name == "ParticleIDs":
            dt = "uint32"
        else:
            dt = "float32"
        if name in self._vector_fields:
            count *= 3
        arr = np.fromfile(f, dtype=dt, count = count)
        if name in self._vector_fields:
            arr = arr.reshape((count/3, 3), order="C")
        return arr.astype("float64")

    def _initialize_index(self, data_file, regions):
        count = sum(data_file.total_particles.values())
        DLE = data_file.ds.domain_left_edge
        DRE = data_file.ds.domain_right_edge
        dx = (DRE - DLE) / 2**_ORDER_MAX
        with open(data_file.filename, "rb") as f:
            # We add on an additionally 4 for the first record.
            f.seek(data_file._position_offset + 4)
            # The first total_particles * 3 values are positions
            pp = np.fromfile(f, dtype = 'float32', count = count*3)
            pp.shape = (count, 3)
        regions.add_data_file(pp, data_file.file_id, data_file.ds.filter_bbox)
        morton = compute_morton(pp[:,0], pp[:,1], pp[:,2], DLE, DRE,
                                data_file.ds.filter_bbox)
        return morton

    def _count_particles(self, data_file):
        npart = dict((self._ptypes[i], v)
            for i, v in enumerate(data_file.header["Npart"]))
        return npart

    # header is 256, but we have 4 at beginning and end for ints
    _field_size = 4
    def _calculate_field_offsets(self, field_list, pcount,
                                 offset, file_size = None):
        # field_list is (ftype, fname) but the blocks are ordered
        # (fname, ftype) in the file.
        pos = offset
        fs = self._field_size
        offsets = {}
        for field in self._fields:
            if not isinstance(field, types.StringTypes):
                field = field[0]
            if not any( (ptype, field) in field_list
                        for ptype in self._ptypes):
                continue
            pos += 4
            any_ptypes = False
            for ptype in self._ptypes:
                if field == "Mass" and ptype not in self.var_mass:
                    continue
                if (ptype, field) not in field_list:
                    continue
                offsets[(ptype, field)] = pos
                any_ptypes = True
                if field in self._vector_fields:
                    pos += 3 * pcount[ptype] * fs
                else:
                    pos += pcount[ptype] * fs
            pos += 4
            if not any_ptypes: pos -= 8
        if file_size is not None:
            if file_size != pos:
                mylog.warning("Your Gadget-2 file may have extra " +
                              "columns or different precision!" +
                              " (%s file vs %s computed)",
                              file_size, pos)
        return offsets

    def _identify_fields(self, domain):
        # We can just look at the particle counts.
        field_list = []
        tp = domain.total_particles
        for i, ptype in enumerate(self._ptypes):
            count = tp[ptype]
            if count == 0: continue
            m = domain.header["Massarr"][i]
            for field in self._fields:
                if isinstance(field, types.TupleType):
                    field, req = field
                    if req is ZeroMass:
                        if m > 0.0 : continue
                    elif req != ptype:
                        continue
                field_list.append((ptype, field))
        return field_list, {}

class IOHandlerTipsyBinary(BaseIOHandler):
    _dataset_type = "tipsy"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")

    _pdtypes = None # dtypes, to be filled in later

    _ptypes = ( "Gas",
                "DarkMatter",
                "Stars" )
    _chunksize = 64*64*64

    _aux_fields = None
    _fields = ( ("Gas", "Mass"),
                ("Gas", "Coordinates"),
                ("Gas", "Velocities"),
                ("Gas", "Density"),
                ("Gas", "Temperature"),
                ("Gas", "Epsilon"),
                ("Gas", "Metals"),
                ("Gas", "Phi"),
                ("DarkMatter", "Mass"),
                ("DarkMatter", "Coordinates"),
                ("DarkMatter", "Velocities"),
                ("DarkMatter", "Epsilon"),
                ("DarkMatter", "Phi"),
                ("Stars", "Mass"),
                ("Stars", "Coordinates"),
                ("Stars", "Velocities"),
                ("Stars", "Metals"),
                ("Stars", "FormationTime"),
                ("Stars", "Epsilon"),
                ("Stars", "Phi")
              )

    def __init__(self, *args, **kwargs):
        self._aux_fields = []
        super(IOHandlerTipsyBinary, self).__init__(*args, **kwargs)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_aux_fields(self, field, mask, data_file):
        """
        Read in auxiliary files from gasoline/pkdgrav.
        This method will automatically detect the format of the file.
        """
        filename = data_file.filename+'.'+field
        dtype = None
        # We need to do some fairly ugly detection to see what format the auxiliary
        # files are in.  They can be either ascii or binary, and the binary files can be
        # either floats, ints, or doubles.  We're going to use a try-catch cascade to
        # determine the format.
        try:#ASCII
            auxdata = np.genfromtxt(filename, skip_header=1)
            if auxdata.size != np.sum(data_file.total_particles.values()):
                print "Error reading auxiliary tipsy file"
                raise RuntimeError
        except ValueError:#binary/xdr
            f = open(filename, 'rb')
            l = struct.unpack(data_file.ds.endian+"i", f.read(4))[0]
            if l != np.sum(data_file.total_particles.values()):
                print "Error reading auxiliary tipsy file"
                raise RuntimeError
            dtype = 'd'
            if field in ('iord', 'igasorder', 'grp'):#These fields are integers
                dtype = 'i'
            try:# If we try loading doubles by default, we can catch an exception and try floats next
                auxdata = np.array(struct.unpack(data_file.ds.endian+(l*dtype), f.read()))
            except struct.error:
                f.seek(4)
                dtype = 'f'
                try:
                    auxdata = np.array(struct.unpack(data_file.ds.endian+(l*dtype), f.read()))
                except struct.error: # None of the binary attempts to read succeeded
                    print "Error reading auxiliary tipsy file"
                    raise RuntimeError

        # Use the mask to slice out the appropriate particle type data
        if mask.size == data_file.total_particles['Gas']:
            return auxdata[:data_file.total_particles['Gas']]
        elif mask.size == data_file.total_particles['DarkMatter']:
            return auxdata[data_file.total_particles['Gas']:-data_file.total_particles['DarkMatter']]
        else:
            return auxdata[-data_file.total_particles['Stars']:]

    def _fill_fields(self, fields, vals, mask, data_file):
        if mask is None:
            size = 0
        else:
            size = mask.sum()
        rv = {}
        for field in fields:
            mylog.debug("Allocating %s values for %s", size, field)
            if field in self._aux_fields: #Read each of the auxiliary fields
                rv[field] = self._read_aux_fields(field, mask, data_file)
            elif field in self._vector_fields:
                rv[field] = np.empty((size, 3), dtype="float64")
                if size == 0: continue
                rv[field][:,0] = vals[field]['x'][mask]
                rv[field][:,1] = vals[field]['y'][mask]
                rv[field][:,2] = vals[field]['z'][mask]
            else:
                rv[field] = np.empty(size, dtype="float64")
                if size == 0: continue
                rv[field][:] = vals[field][mask]
            if field == "Coordinates":
                eps = np.finfo(rv[field].dtype).eps
                for i in range(3):
                  rv[field][:,i] = np.clip(rv[field][:,i],
                      self.domain_left_edge[i] + eps,
                      self.domain_right_edge[i] - eps)
        return rv


    def _read_particle_coords(self, chunks, ptf):
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items(), key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                total = 0
                while total < tp[ptype]:
                    p = np.fromfile(f, self._pdtypes[ptype],
                            count=min(self._chunksize, tp[ptype] - total))
                    total += p.size
                    d = [p["Coordinates"][ax].astype("float64") for ax in 'xyz']
                    del p
                    yield ptype, d

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            poff = data_file.field_offsets
            tp = data_file.total_particles
            f = open(data_file.filename, "rb")
            for ptype, field_list in sorted(ptf.items(), key=lambda a: poff[a[0]]):
                f.seek(poff[ptype], os.SEEK_SET)
                total = 0
                while total < tp[ptype]:
                    p = np.fromfile(f, self._pdtypes[ptype],
                        count=min(self._chunksize, tp[ptype] - total))
                    total += p.size
                    mask = selector.select_points(
                        p["Coordinates"]['x'].astype("float64"),
                        p["Coordinates"]['y'].astype("float64"),
                        p["Coordinates"]['z'].astype("float64"), 0.0)
                    if mask is None: continue
                    tf = self._fill_fields(field_list, p, mask, data_file)
                    for field in field_list:
                        yield (ptype, field), tf.pop(field)
            f.close()

    def _update_domain(self, data_file):
        '''
        This method is used to determine the size needed for a box that will
        bound the particles.  It simply finds the largest position of the
        whole set of particles, and sets the domain to +/- that value.
        '''
        ds = data_file.ds
        ind = 0
        # Check to make sure that the domain hasn't already been set
        # by the parameter file
        if np.all(np.isfinite(ds.domain_left_edge)) and np.all(np.isfinite(ds.domain_right_edge)):
            return
        with open(data_file.filename, "rb") as f:
            ds.domain_left_edge = 0
            ds.domain_right_edge = 0
            f.seek(ds._header_offset)
            mi =   np.array([1e30, 1e30, 1e30], dtype="float64")
            ma =  -np.array([1e30, 1e30, 1e30], dtype="float64")
            for iptype, ptype in enumerate(self._ptypes):
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0: continue
                start, stop = ind, ind + count
                while ind < stop:
                    c = min(CHUNKSIZE, stop - ind)
                    pp = np.fromfile(f, dtype = self._pdtypes[ptype],
                                     count = c)
                    eps = np.finfo(pp["Coordinates"]["x"].dtype).eps
                    np.minimum(mi, [pp["Coordinates"]["x"].min(),
                                    pp["Coordinates"]["y"].min(),
                                    pp["Coordinates"]["z"].min()], mi)
                    np.maximum(ma, [pp["Coordinates"]["x"].max(),
                                    pp["Coordinates"]["y"].max(),
                                    pp["Coordinates"]["z"].max()], ma)
                    ind += c
        # We extend by 1%.
        DW = ma - mi
        mi -= 0.01 * DW
        ma += 0.01 * DW
        ds.domain_left_edge = ds.arr(mi, 'code_length')
        ds.domain_right_edge = ds.arr(ma, 'code_length')
        ds.domain_width = DW = ds.domain_right_edge - ds.domain_left_edge
        ds.unit_registry.add("unitary", float(DW.max() * DW.units.cgs_value),
                                 DW.units.dimensions)

    def _initialize_index(self, data_file, regions):
        ds = data_file.ds
        morton = np.empty(sum(data_file.total_particles.values()),
                          dtype="uint64")
        ind = 0
        DLE, DRE = ds.domain_left_edge, ds.domain_right_edge
        dx = (DRE - DLE) / (2**_ORDER_MAX)
        self.domain_left_edge = DLE.in_units("code_length").ndarray_view()
        self.domain_right_edge = DRE.in_units("code_length").ndarray_view()
        with open(data_file.filename, "rb") as f:
            f.seek(ds._header_offset)
            for iptype, ptype in enumerate(self._ptypes):
                # We'll just add the individual types separately
                count = data_file.total_particles[ptype]
                if count == 0: continue
                start, stop = ind, ind + count
                while ind < stop:
                    c = min(CHUNKSIZE, stop - ind)
                    pp = np.fromfile(f, dtype = self._pdtypes[ptype],
                                     count = c)
                    mis = np.empty(3, dtype="float64")
                    mas = np.empty(3, dtype="float64")
                    for axi, ax in enumerate('xyz'):
                        mi = pp["Coordinates"][ax].min()
                        ma = pp["Coordinates"][ax].max()
                        mylog.debug("Spanning: %0.3e .. %0.3e in %s", mi, ma, ax)
                        mis[axi] = mi
                        mas[axi] = ma
                    pos = np.empty((pp.size, 3), dtype="float64")
                    for i, ax in enumerate("xyz"):
                        eps = np.finfo(pp["Coordinates"][ax].dtype).eps
                        pos[:,i] = pp["Coordinates"][ax]
                    regions.add_data_file(pos, data_file.file_id,
                                          data_file.ds.filter_bbox)
                    morton[ind:ind+c] = compute_morton(
                        pos[:,0], pos[:,1], pos[:,2],
                        DLE, DRE, data_file.ds.filter_bbox)
                    ind += c
        mylog.info("Adding %0.3e particles", morton.size)
        return morton

    def _count_particles(self, data_file):
        npart = {
            "Gas": data_file.ds.parameters['nsph'],
            "Stars": data_file.ds.parameters['nstar'],
            "DarkMatter": data_file.ds.parameters['ndark']
        }
        return npart

    @classmethod
    def _compute_dtypes(cls, field_dtypes, endian = "<"):
        pds = {}
        for ptype, field in cls._fields:
            dtbase = field_dtypes.get(field, 'f')
            ff = "%s%s" % (endian, dtbase)
            if field in cls._vector_fields:
                dt = (field, [('x', ff), ('y', ff), ('z', ff)])
            else:
                dt = (field, ff)
            pds.setdefault(ptype, []).append(dt)
        pdtypes = {}
        for ptype in pds:
            pdtypes[ptype] = np.dtype(pds[ptype])
        return pdtypes

    def _create_dtypes(self, data_file):
        # We can just look at the particle counts.
        self._header_offset = data_file.ds._header_offset
        self._pdtypes = {}
        pds = {}
        field_list = []
        tp = data_file.total_particles
        aux_filenames = glob.glob(data_file.filename+'.*') # Find out which auxiliaries we have
        self._aux_fields = [f[1+len(data_file.filename):] for f in aux_filenames]
        self._pdtypes = self._compute_dtypes(data_file.ds._field_dtypes,
                                             data_file.ds.endian)
        for ptype, field in self._fields:
            if tp[ptype] == 0:
                # We do not want out _pdtypes to have empty particles.
                self._pdtypes.pop(ptype, None)
                continue
            field_list.append((ptype, field))
        if any(["Gas"==f[0] for f in field_list]): #Add the auxiliary fields to each ptype we have
            field_list += [("Gas",a) for a in self._aux_fields]
        if any(["DarkMatter"==f[0] for f in field_list]):
            field_list += [("DarkMatter",a) for a in self._aux_fields]
        if any(["Stars"==f[0] for f in field_list]):
            field_list += [("Stars",a) for a in self._aux_fields]
        self._field_list = field_list
        return self._field_list

    def _identify_fields(self, data_file):
        return self._field_list, {}

    def _calculate_particle_offsets(self, data_file):
        field_offsets = {}
        pos = data_file.ds._header_offset
        for ptype in self._ptypes:
            field_offsets[ptype] = pos
            if data_file.total_particles[ptype] == 0: continue
            size = self._pdtypes[ptype].itemsize
            pos += data_file.total_particles[ptype] * size
        return field_offsets

class IOHandlerHTTPStream(BaseIOHandler):
    _dataset_type = "http_particle_stream"
    _vector_fields = ("Coordinates", "Velocity", "Velocities")

    def __init__(self, ds):
        if requests is None:
            raise RuntimeError
        self._url = ds.base_url
        # This should eventually manage the IO and cache it
        self.total_bytes = 0
        super(IOHandlerHTTPStream, self).__init__(ds)

    def _open_stream(self, data_file, field):
        # This does not actually stream yet!
        ftype, fname = field
        s = "%s/%s/%s/%s" % (self._url,
            data_file.file_id, ftype, fname)
        mylog.info("Loading URL %s", s)
        resp = requests.get(s)
        if resp.status_code != 200:
            raise RuntimeError
        self.total_bytes += len(resp.content)
        return resp.content

    def _identify_fields(self, data_file):
        f = []
        for ftype, fname in self.ds.parameters["field_list"]:
            f.append((str(ftype), str(fname)))
        return f, {}

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            for ptype in ptf:
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0]/3.0, 3)
                yield ptype, (c[:,0], c[:,1], c[:,2])

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            for ptype, field_list in sorted(ptf.items()):
                s = self._open_stream(data_file, (ptype, "Coordinates"))
                c = np.frombuffer(s, dtype="float64")
                c.shape = (c.shape[0]/3.0, 3)
                mask = selector.select_points(
                            c[:,0], c[:,1], c[:,2], 0.0)
                del c
                if mask is None: continue
                for field in field_list:
                    s = self._open_stream(data_file, (ptype, field))
                    c = np.frombuffer(s, dtype="float64")
                    if field in self._vector_fields:
                        c.shape = (c.shape[0]/3.0, 3)
                    data = c[mask, ...]
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        header = self.ds.parameters
        ptypes = header["particle_count"][data_file.file_id].keys()
        pcount = sum(header["particle_count"][data_file.file_id].values())
        morton = np.empty(pcount, dtype='uint64')
        ind = 0
        for ptype in ptypes:
            s = self._open_stream(data_file, (ptype, "Coordinates"))
            c = np.frombuffer(s, dtype="float64")
            c.shape = (c.shape[0]/3.0, 3)
            regions.add_data_file(c, data_file.file_id,
                                  data_file.ds.filter_bbox)
            morton[ind:ind+c.shape[0]] = compute_morton(
                c[:,0], c[:,1], c[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge,
                data_file.ds.filter_bbox)
            ind += c.shape[0]
        return morton

    def _count_particles(self, data_file):
        return self.ds.parameters["particle_count"][data_file.file_id]
