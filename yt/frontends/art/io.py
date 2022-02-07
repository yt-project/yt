import os
import os.path
from collections import defaultdict
from functools import partial

import numpy as np

from yt.frontends.art.definitions import (
    hydro_struct,
    particle_fields,
    particle_star_fields,
    star_struct,
)
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.fortran_utils import read_vector, skip
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


class IOHandlerART(BaseIOHandler):
    _dataset_type = "art"
    tb, ages = None, None
    cache = None
    masks = None
    caching = True

    def __init__(self, *args, **kwargs):
        self.cache = {}
        self.masks = {}
        super().__init__(*args, **kwargs)
        self.ws = self.ds.parameters["wspecies"]
        self.ls = self.ds.parameters["lspecies"]
        self.file_particle = self.ds._file_particle_data
        self.file_stars = self.ds._file_particle_stars
        self.Nrow = self.ds.parameters["Nrow"]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Chunks in this case will have affiliated domain subset objects
        # Each domain subset will contain a hydro_offset array, which gives
        # pointers to level-by-level hydro information
        tr = defaultdict(list)
        cp = 0
        for chunk in chunks:
            for subset in chunk.objs:
                # Now we read the entire thing
                f = open(subset.domain.ds._file_amr, "rb")
                # This contains the boundary information, so we skim through
                # and pick off the right vectors
                rv = subset.fill(f, fields, selector)
                for ft, f in fields:
                    d = rv.pop(f)
                    mylog.debug(
                        "Filling %s with %s (%0.3e %0.3e) (%s:%s)",
                        f,
                        d.size,
                        d.min(),
                        d.max(),
                        cp,
                        cp + d.size,
                    )
                    tr[(ft, f)].append(d)
                cp += d.size
        d = {}
        for field in fields:
            d[field] = np.concatenate(tr.pop(field))
        return d

    def _get_mask(self, selector, ftype):
        key = (selector, ftype)
        if key in self.masks.keys() and self.caching:
            return self.masks[key]
        pstr = "particle_position_%s"
        x, y, z = (self._get_field((ftype, pstr % ax)) for ax in "xyz")
        mask = selector.select_points(x, y, z, 0.0)
        if self.caching:
            self.masks[key] = mask
            return self.masks[key]
        else:
            return mask

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for _chunk in chunks:
            for ptype in sorted(ptf):
                x = self._get_field((ptype, "particle_position_x"))
                y = self._get_field((ptype, "particle_position_y"))
                z = self._get_field((ptype, "particle_position_z"))
                yield ptype, (x, y, z), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for _chunk in chunks:
            for ptype, field_list in sorted(ptf.items()):
                x = self._get_field((ptype, "particle_position_x"))
                y = self._get_field((ptype, "particle_position_y"))
                z = self._get_field((ptype, "particle_position_z"))
                mask = selector.select_points(x, y, z, 0.0)
                if mask is None:
                    continue
                for field in field_list:
                    data = self._get_field((ptype, field))
                    yield (ptype, field), data[mask]

    def _get_field(self, field):
        if field in self.cache.keys() and self.caching:
            mylog.debug("Cached %s", str(field))
            return self.cache[field]
        mylog.debug("Reading %s", str(field))
        tr = {}
        ftype, fname = field
        ptmax = self.ws[-1]
        pbool, idxa, idxb = _determine_field_size(self.ds, ftype, self.ls, ptmax)
        npa = idxb - idxa
        sizes = np.diff(np.concatenate(([0], self.ls)))
        rp = partial(
            read_particles, self.file_particle, self.Nrow, idxa=idxa, idxb=idxb
        )
        for ax in "xyz":
            if fname.startswith(f"particle_position_{ax}"):
                dd = self.ds.domain_dimensions[0]
                off = 1.0 / dd
                tr[field] = rp(fields=[ax])[0] / dd - off
            if fname.startswith(f"particle_velocity_{ax}"):
                (tr[field],) = rp(fields=["v" + ax])
        if fname.startswith("particle_mass"):
            a = 0
            data = np.zeros(npa, dtype="f8")
            for ptb, size, m in zip(pbool, sizes, self.ws):
                if ptb:
                    data[a : a + size] = m
                    a += size
            tr[field] = data
        elif fname == "particle_index":
            tr[field] = np.arange(idxa, idxb)
        elif fname == "particle_type":
            a = 0
            data = np.zeros(npa, dtype="int64")
            for i, (ptb, size) in enumerate(zip(pbool, sizes)):
                if ptb:
                    data[a : a + size] = i
                    a += size
            tr[field] = data
        if pbool[-1] and fname in particle_star_fields:
            data = read_star_field(self.file_stars, field=fname)
            temp = tr.get(field, np.zeros(npa, "f8"))
            nstars = self.ls[-1] - self.ls[-2]
            if nstars > 0:
                temp[-nstars:] = data
            tr[field] = temp
        if fname == "particle_creation_time":
            self.tb, self.ages, data = interpolate_ages(
                tr[field][-nstars:],
                self.file_stars,
                self.tb,
                self.ages,
                self.ds.current_time,
            )
            temp = tr.get(field, np.zeros(npa, "f8"))
            temp[-nstars:] = data
            tr[field] = temp
            del data
        # We check again, after it's been filled
        if fname.startswith("particle_mass"):
            # We now divide by NGrid in order to make this match up.  Note that
            # this means that even when requested in *code units*, we are
            # giving them as modified by the ng value.  This only works for
            # dark_matter -- stars are regular matter.
            tr[field] /= self.ds.domain_dimensions.prod()
        if tr == {}:
            tr = {f: np.array([]) for f in [field]}
        if self.caching:
            self.cache[field] = tr[field]
            return self.cache[field]
        else:
            return tr[field]


class IOHandlerDarkMatterART(IOHandlerART):
    _dataset_type = "dm_art"

    def _count_particles(self, data_file):
        return {
            k: self.ds.parameters["lspecies"][i]
            for i, k in enumerate(self.ds.particle_types_raw)
        }

    def _identify_fields(self, domain):
        field_list = []
        self.particle_field_list = [f for f in particle_fields]
        for ptype in self.ds.particle_types_raw:
            for pfield in self.particle_field_list:
                pfn = (ptype, pfield)
                field_list.append(pfn)
        return field_list, {}

    def _get_field(self, field):
        if field in self.cache.keys() and self.caching:
            mylog.debug("Cached %s", str(field))
            return self.cache[field]
        mylog.debug("Reading %s", str(field))
        tr = {}
        ftype, fname = field
        ptmax = self.ws[-1]
        pbool, idxa, idxb = _determine_field_size(self.ds, ftype, self.ls, ptmax)
        npa = idxb - idxa
        sizes = np.diff(np.concatenate(([0], self.ls)))
        rp = partial(
            read_particles, self.file_particle, self.Nrow, idxa=idxa, idxb=idxb
        )
        for ax in "xyz":
            if fname.startswith(f"particle_position_{ax}"):
                # This is not the same as domain_dimensions
                dd = self.ds.parameters["ng"]
                off = 1.0 / dd
                tr[field] = rp(fields=[ax])[0] / dd - off
            if fname.startswith(f"particle_velocity_{ax}"):
                (tr[field],) = rp(["v" + ax])
        if fname.startswith("particle_mass"):
            a = 0
            data = np.zeros(npa, dtype="f8")
            for ptb, size, m in zip(pbool, sizes, self.ws):
                if ptb:
                    data[a : a + size] = m
                    a += size
            tr[field] = data
        elif fname == "particle_index":
            tr[field] = np.arange(idxa, idxb)
        elif fname == "particle_type":
            a = 0
            data = np.zeros(npa, dtype="int64")
            for i, (ptb, size) in enumerate(zip(pbool, sizes)):
                if ptb:
                    data[a : a + size] = i
                    a += size
            tr[field] = data
        # We check again, after it's been filled
        if fname.startswith("particle_mass"):
            # We now divide by NGrid in order to make this match up.  Note that
            # this means that even when requested in *code units*, we are
            # giving them as modified by the ng value.  This only works for
            # dark_matter -- stars are regular matter.
            tr[field] /= self.ds.domain_dimensions.prod()
        if tr == {}:
            tr[field] = np.array([])
        if self.caching:
            self.cache[field] = tr[field]
            return self.cache[field]
        else:
            return tr[field]

    def _yield_coordinates(self, data_file):
        for ptype in self.ds.particle_types_raw:
            x = self._get_field((ptype, "particle_position_x"))
            y = self._get_field((ptype, "particle_position_y"))
            z = self._get_field((ptype, "particle_position_z"))

            yield ptype, np.stack((x, y, z), axis=-1)


def _determine_field_size(pf, field, lspecies, ptmax):
    pbool = np.zeros(len(lspecies), dtype="bool")
    idxas = np.concatenate(
        (
            [
                0,
            ],
            lspecies[:-1],
        )
    )
    idxbs = lspecies
    if "specie" in field:
        index = int(field.replace("specie", ""))
        pbool[index] = True
    else:
        raise RuntimeError
    idxa, idxb = idxas[pbool][0], idxbs[pbool][-1]
    return pbool, idxa, idxb


def interpolate_ages(
    data, file_stars, interp_tb=None, interp_ages=None, current_time=None
):
    if interp_tb is None:
        t_stars, a_stars = read_star_field(file_stars, field="t_stars")
        # timestamp of file should match amr timestamp
        if current_time:
            tdiff = YTQuantity(b2t(t_stars), "Gyr") - current_time.in_units("Gyr")
            if np.abs(tdiff) > 1e-4:
                mylog.info("Timestamp mismatch in star particle header: %s", tdiff)
        mylog.info("Interpolating ages")
        interp_tb, interp_ages = b2t(data)
        interp_tb = YTArray(interp_tb, "Gyr")
        interp_ages = YTArray(interp_ages, "Gyr")
    temp = np.interp(data, interp_tb, interp_ages)
    return interp_tb, interp_ages, temp


def _read_art_level_info(
    f, level_oct_offsets, level, coarse_grid=128, ncell0=None, root_level=None
):
    pos = f.tell()
    f.seek(level_oct_offsets[level])
    # Get the info for this level, skip the rest
    junk, nLevel, iOct = read_vector(f, "i", ">")

    # fortran indices start at 1

    # Skip all the oct index data
    le = np.zeros((nLevel, 3), dtype="int64")
    fl = np.ones((nLevel, 6), dtype="int64")
    iocts = np.zeros(nLevel + 1, dtype="int64")
    idxa, idxb = 0, 0
    chunk = int(1e6)  # this is ~111MB for 15 dimensional 64 bit arrays
    left = nLevel
    while left > 0:
        this_chunk = min(chunk, left)
        idxb = idxa + this_chunk
        data = np.fromfile(f, dtype=">i", count=this_chunk * 15)
        data = data.reshape(this_chunk, 15)
        left -= this_chunk
        le[idxa:idxb, :] = data[:, 1:4]
        fl[idxa:idxb, 1] = np.arange(idxa, idxb)
        # pad byte is last, LL2, then ioct right before it
        iocts[idxa:idxb] = data[:, -3]
        idxa = idxa + this_chunk
    del data

    # emulate fortran code
    #     do ic1 = 1 , nLevel
    #       read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
    # &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct),
    # &                iOctLL2(iOct)
    #       iOct = iOctLL1(iOct)

    # ioct always represents the index of the next variable
    # not the current, so shift forward one index
    # the last index isn't used
    iocts[1:] = iocts[:-1]  # shift
    iocts = iocts[:nLevel]  # chop off the last, unused, index
    iocts[0] = iOct  # starting value

    # now correct iocts for fortran indices start @ 1
    iocts = iocts - 1

    assert np.unique(iocts).shape[0] == nLevel

    # left edges are expressed as if they were on
    # level 15, so no matter what level max(le)=2**15
    # correct to the yt convention
    # le = le/2**(root_level-1-level)-1

    # try to find the root_level first
    def cfc(root_level, level, le):
        d_x = 1.0 / (2.0 ** (root_level - level + 1))
        fc = (d_x * le) - 2 ** (level - 1)
        return fc

    if root_level is None:
        root_level = np.floor(np.log2(le.max() * 1.0 / coarse_grid))
        root_level = root_level.astype("int64")
        for _ in range(10):
            fc = cfc(root_level, level, le)
            go = np.diff(np.unique(fc)).min() < 1.1
            if go:
                break
            root_level += 1
    else:
        fc = cfc(root_level, level, le)
    unitary_center = fc / (coarse_grid * 2.0 ** (level - 1))
    assert np.all(unitary_center < 1.0)

    # again emulate the fortran code
    # This is all for calculating child oct locations
    # iC_ = iC + nbshift
    # iO = ishft ( iC_ , - ndim )
    # id = ishft ( 1, MaxLevel - iOctLv(iO) )
    # j  = iC_ + 1 - ishft( iO , ndim )
    # Posx   = d_x * (iOctPs(1,iO) + sign ( id , idelta(j,1) ))
    # Posy   = d_x * (iOctPs(2,iO) + sign ( id , idelta(j,2) ))
    # Posz   = d_x * (iOctPs(3,iO) + sign ( id , idelta(j,3) ))
    # idelta = [[-1,  1, -1,  1, -1,  1, -1,  1],
    #           [-1, -1,  1,  1, -1, -1,  1,  1],
    #           [-1, -1, -1, -1,  1,  1,  1,  1]]
    # idelta = np.array(idelta)
    # if ncell0 is None:
    #     ncell0 = coarse_grid**3
    # nchild = 8
    # ndim = 3
    # nshift = nchild -1
    # nbshift = nshift - ncell0
    # iC = iocts #+ nbshift
    # iO = iC >> ndim #possibly >>
    # id = 1 << (root_level - level)
    # j = iC + 1 - ( iO << 3)
    # delta = np.abs(id)*idelta[:,j-1]

    # try without the -1
    # le = le/2**(root_level+1-level)
    # now read the hvars and vars arrays
    # we are looking for iOctCh
    # we record if iOctCh is >0, in which it is subdivided
    # iOctCh  = np.zeros((nLevel+1,8),dtype='bool')
    f.seek(pos)
    return unitary_center, fl, iocts, nLevel, root_level


def get_ranges(
    skip, count, field, words=6, real_size=4, np_per_page=4096**2, num_pages=1
):
    # translate every particle index into a file position ranges
    ranges = []
    arr_size = np_per_page * real_size
    idxa, idxb = 0, 0
    posa, posb = 0, 0
    for _page in range(num_pages):
        idxb += np_per_page
        for i, fname in enumerate(["x", "y", "z", "vx", "vy", "vz"]):
            posb += arr_size
            if i == field or fname == field:
                if skip < np_per_page and count > 0:
                    left_in_page = np_per_page - skip
                    this_count = min(left_in_page, count)
                    count -= this_count
                    start = posa + skip * real_size
                    end = posa + this_count * real_size
                    ranges.append((start, this_count))
                    skip = 0
                    assert end <= posb
                else:
                    skip -= np_per_page
            posa += arr_size
        idxa += np_per_page
    assert count == 0
    return ranges


def read_particles(file, Nrow, idxa, idxb, fields):
    words = 6  # words (reals) per particle: x,y,z,vx,vy,vz
    real_size = 4  # for file_particle_data; not always true?
    np_per_page = Nrow**2  # defined in ART a_setup.h, # of particles/page
    num_pages = os.path.getsize(file) // (real_size * words * np_per_page)
    fh = open(file)
    skip, count = idxa, idxb - idxa
    kwargs = dict(
        words=words, real_size=real_size, np_per_page=np_per_page, num_pages=num_pages
    )
    arrs = []
    for field in fields:
        ranges = get_ranges(skip, count, field, **kwargs)
        data = None
        for seek, this_count in ranges:
            fh.seek(seek)
            temp = np.fromfile(fh, count=this_count, dtype=">f4")
            if data is None:
                data = temp
            else:
                data = np.concatenate((data, temp))
        arrs.append(data.astype("f8"))
    fh.close()
    return arrs


def read_star_field(file, field=None):
    data = {}
    with open(file, "rb") as fh:
        for dtype, variables in star_struct:
            found = (
                isinstance(variables, tuple) and field in variables
            ) or field == variables
            if found:
                data[field] = read_vector(fh, dtype[1], dtype[0])
            else:
                skip(fh, endian=">")
    return data.pop(field)


def _read_child_mask_level(f, level_child_offsets, level, nLevel, nhydro_vars):
    f.seek(level_child_offsets[level])
    ioctch = np.zeros(nLevel, dtype="uint8")
    idc = np.zeros(nLevel, dtype="int32")

    chunk = int(1e6)
    left = nLevel
    width = nhydro_vars + 6
    a, b = 0, 0
    while left > 0:
        chunk = min(chunk, left)
        b += chunk
        arr = np.fromfile(f, dtype=">i", count=chunk * width)
        arr = arr.reshape((width, chunk), order="F")
        assert np.all(arr[0, :] == arr[-1, :])  # pads must be equal
        idc[a:b] = arr[1, :] - 1  # fix fortran indexing
        ioctch[a:b] = arr[2, :] == 0  # if it is above zero, then refined available
        # zero in the mask means there is refinement available
        a = b
        left -= chunk
    assert left == 0
    return idc, ioctch


nchem = 8 + 2
dtyp = np.dtype(f">i4,>i8,>i8,>{nchem}f4,>2f4,>i4")


def _read_child_level(
    f,
    level_child_offsets,
    level_oct_offsets,
    level_info,
    level,
    fields,
    domain_dimensions,
    ncell0,
    nhydro_vars=10,
    nchild=8,
    noct_range=None,
):
    # emulate the fortran code for reading cell data
    # read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar),
    #    &                 (var(i,idc), i=2,3)
    # contiguous 8-cell sections are for the same oct;
    # ie, we don't write out just the 0 cells, then the 1 cells
    # optionally, we only read noct_range to save memory
    left_index, fl, octs, nocts, root_level = _read_art_level_info(
        f, level_oct_offsets, level, coarse_grid=domain_dimensions[0]
    )
    if noct_range is None:
        nocts = level_info[level]
        ncells = nocts * 8
        f.seek(level_child_offsets[level])
        arr = np.fromfile(f, dtype=hydro_struct, count=ncells)
        assert np.all(arr["pad1"] == arr["pad2"])  # pads must be equal
        # idc = np.argsort(arr['idc']) #correct fortran indices
        # translate idc into icell, and then to iOct
        icell = (arr["idc"] >> 3) << 3
        iocts = (icell - ncell0) / nchild  # without a F correction, there's a +1
        # assert that the children are read in the same order as the octs
        assert np.all(octs == iocts[::nchild])
    else:
        start, end = noct_range
        nocts = min(end - start, level_info[level])
        end = start + nocts
        ncells = nocts * 8
        skip = np.dtype(hydro_struct).itemsize * start * 8
        f.seek(level_child_offsets[level] + skip)
        arr = np.fromfile(f, dtype=hydro_struct, count=ncells)
        assert np.all(arr["pad1"] == arr["pad2"])  # pads must be equal
    source = {}
    for field in fields:
        sh = (nocts, 8)
        source[field] = np.reshape(arr[field], sh, order="C").astype("float64")
    return source


def _read_root_level(f, level_offsets, level_info, nhydro_vars=10):
    nocts = level_info[0]
    f.seek(level_offsets[0])  # Ditch the header
    hvar = read_vector(f, "f", ">")
    var = read_vector(f, "f", ">")
    hvar = hvar.reshape((nhydro_vars, nocts * 8), order="F")
    var = var.reshape((2, nocts * 8), order="F")
    arr = np.concatenate((hvar, var))
    return arr


# All of these functions are to convert from hydro time var to
# proper time
sqrt = np.sqrt
sign = np.sign


def find_root(f, a, b, tol=1e-6):
    c = (a + b) / 2.0
    last = -np.inf
    assert sign(f(a)) != sign(f(b))
    while np.abs(f(c) - last) > tol:
        last = f(c)
        if sign(last) == sign(f(b)):
            b = c
        else:
            a = c
        c = (a + b) / 2.0
    return c


def quad(fintegrand, xmin, xmax, n=1e4):
    spacings = np.logspace(np.log10(xmin), np.log10(xmax), num=int(n))
    integrand_arr = fintegrand(spacings)
    val = np.trapz(integrand_arr, dx=np.diff(spacings))
    return val


def a2b(at, Om0=0.27, Oml0=0.73, h=0.700):
    def f_a2b(x):
        val = 0.5 * sqrt(Om0) / x**3.0
        val /= sqrt(Om0 / x**3.0 + Oml0 + (1.0 - Om0 - Oml0) / x**2.0)
        return val

    # val, err = si.quad(f_a2b,1,at)
    val = quad(f_a2b, 1, at)
    return val


def b2a(bt, **kwargs):
    # converts code time into expansion factor
    # if Om0 ==1and OmL == 0 then b2a is (1 / (1-td))**2
    # if bt < -190.0 or bt > -.10:  raise 'bt outside of range'
    def f_b2a(at):
        return a2b(at, **kwargs) - bt

    return find_root(f_b2a, 1e-4, 1.1)
    # return so.brenth(f_b2a,1e-4,1.1)
    # return brent.brent(f_b2a)


def a2t(at, Om0=0.27, Oml0=0.73, h=0.700):
    def integrand(x):
        return 1.0 / (x * sqrt(Oml0 + Om0 * x**-3.0))

    # current_time,err = si.quad(integrand,0.0,at,epsabs=1e-6,epsrel=1e-6)
    current_time = quad(integrand, 1e-4, at)
    # spacings = np.logspace(-5,np.log10(at),num=int(1e5))
    # integrand_arr = integrand(spacings)
    # current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
    current_time *= 9.779 / h
    return current_time


def b2t(tb, n=1e2, logger=None, **kwargs):
    tb = np.array(tb)
    if isinstance(tb, float):
        return a2t(b2a(tb))
    if tb.shape == ():
        return a2t(b2a(tb))
    if len(tb) < n:
        n = len(tb)
    tbs = -1.0 * np.logspace(np.log10(-tb.min()), np.log10(-tb.max()), n)
    ages = []
    for i, tbi in enumerate(tbs):
        ages += (a2t(b2a(tbi)),)
        if logger:
            logger(i)
    ages = np.array(ages)
    return tbs, ages
