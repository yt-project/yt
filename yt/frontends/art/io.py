"""
ART-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Chris Moody <matthewturk@gmail.com>
Affiliation: UCSC
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import struct
import os
import os.path

from yt.utilities.io_handler import \
    BaseIOHandler
import yt.utilities.lib as au
from yt.utilities.fortran_utils import *
from yt.utilities.logger import ytLogger as mylog
from yt.frontends.art.definitions import *
from yt.utilities.physical_constants import sec_per_year


class IOHandlerART(BaseIOHandler):
    _data_style = "art"
    tb, ages = None, None

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Chunks in this case will have affiliated domain subset objects
        # Each domain subset will contain a hydro_offset array, which gives
        # pointers to level-by-level hydro information
        tr = dict((f, np.empty(size, dtype='float64')) for f in fields)
        cp = 0
        for chunk in chunks:
            for subset in chunk.objs:
                # Now we read the entire thing
                f = open(subset.domain.pf.file_amr, "rb")
                # This contains the boundary information, so we skim through
                # and pick off the right vectors
                rv = subset.fill(f, fields)
                for ft, f in fields:
                    mylog.debug("Filling %s with %s (%0.3e %0.3e) (%s:%s)",
                                f, subset.cell_count, rv[f].min(), rv[f].max(),
                                cp, cp+subset.cell_count)
                    tr[(ft, f)][cp:cp+subset.cell_count] = rv.pop(f)
                cp += subset.cell_count
        return tr

    def _read_particle_selection(self, chunks, selector, fields):
        # ignore chunking; we have no particle chunk system
        chunk = chunks.next()
        level = chunk.objs[0].domain.domain_level 
        # only chunk out particles on level zero
        if level > 0: 
            tr = {}
            for field in fields:
                tr[field] = np.array([],'f8')
            return tr
        pf = chunk.objs[0].domain.pf
        masks = {}
        ws, ls = pf.parameters["wspecies"], pf.parameters["lspecies"]
        sizes = np.diff(np.concatenate(([0], ls)))
        ptmax = ws[-1]
        npt = ls[-1]
        nstars = ls[-1]-ls[-2]
        file_particle = pf.file_particle_data
        file_stars = pf.file_particle_stars
        tr = {}
        ftype_old = None
        for field in fields:
            ftype, fname = field
            pbool, idxa, idxb = _determine_field_size(pf, ftype, ls, ptmax)
            npa = idxb-idxa
            if not ftype_old == ftype:
                pos, vel = read_particles(file_particle, pf.parameters['Nrow'],
                                          dd=pf.domain_dimensions,
                                          idxa=idxa, idxb=idxb)
                pos, vel = pos.astype('float64'), vel.astype('float64')
                pos -= 1.0/pf.domain_dimensions[0]
                mask = selector.select_points(pos[:, 0], pos[:, 1], pos[:, 2])
                size = mask.sum()
            for i, ax in enumerate('xyz'):
                if fname.startswith("particle_position_%s" % ax):
                    tr[field] = pos[:, i]
                if fname.startswith("particle_velocity_%s" % ax):
                    tr[field] = vel[:, i]
            if fname == "particle_mass":
                a = 0
                data = np.zeros(npa, dtype='f8')
                for ptb, size, m in zip(pbool, sizes, ws):
                    if ptb:
                        data[a:a+size] = m
                        a += size
                tr[field] = data
            elif fname == "particle_index":
                tr[field] = np.arange(idxa, idxb).astype('int64')
            elif fname == "particle_type":
                a = 0
                data = np.zeros(npa, dtype='int')
                for i, (ptb, size) in enumerate(zip(pbool, sizes)):
                    if ptb:
                        data[a:a+size] = i
                        a += size
                tr[field] = data
            if pbool[-1] and fname in particle_star_fields:
                data = read_star_field(file_stars, field=fname)
                temp = tr.get(field, np.zeros(npa, 'f8'))
                if nstars > 0:
                    temp[-nstars:] = data
                tr[field] = temp
            if fname == "particle_creation_time":
                self.tb, self.ages, data = interpolate_ages(
                    tr[field][-nstars:],
                    file_stars,
                    self.tb,
                    self.ages,
                    pf.current_time)
                temp = tr.get(field, np.zeros(npa, 'f8'))
                temp[-nstars:] = data
                tr[field]=temp
                del data
            tr[field] = tr[field][mask]
            ftype_old = ftype
        return tr


def _determine_field_size(pf, field, lspecies, ptmax):
    pbool = np.zeros(len(lspecies), dtype="bool")
    idxas = np.concatenate(([0, ], lspecies[:-1]))
    idxbs = lspecies
    if "specie" in field:
        index = int(field.replace("specie", ""))
        pbool[index] = True
    elif field == "stars":
        pbool[-1] = True
    elif field == "darkmatter":
        pbool[0:-1] = True
    else:
        pbool[:] = True
    idxa, idxb = idxas[pbool][0], idxbs[pbool][-1]
    return pbool, idxa, idxb


def interpolate_ages(data, file_stars, interp_tb=None, interp_ages=None,
                     current_time=None):
    if interp_tb is None:
        tdum, adum = read_star_field(file_stars,
                                     field="tdum")
        # timestamp of file should match amr timestamp
        if current_time:
            tdiff = b2t(tdum)-current_time/(sec_per_year*1e9)
            if np.abs(tdiff) < 1e-4:
                mylog.info("Timestamp mismatch in star " +
                           "particle header")
        mylog.info("Interpolating ages")
        interp_tb, interp_ages = b2t(data)
    temp = np.interp(data, interp_tb, interp_ages)
    temp *= 1.0e9*sec_per_year
    return interp_tb, interp_ages, temp


def _count_art_octs(f, offset,
                    MinLev, MaxLevelNow):
    level_oct_offsets = [0, ]
    level_child_offsets = [0, ]
    f.seek(offset)
    nchild, ntot = 8, 0
    Level = np.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    iNOLL = np.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    iHOLL = np.zeros(MaxLevelNow+1 - MinLev, dtype='int64')
    for Lev in xrange(MinLev + 1, MaxLevelNow+1):
        level_oct_offsets.append(f.tell())

        # Get the info for this level, skip the rest
        # print "Reading oct tree data for level", Lev
        # print 'offset:',f.tell()
        Level[Lev], iNOLL[Lev], iHOLL[Lev] = read_vector(f, 'i', '>')
        # print 'Level %i : '%Lev, iNOLL
        # print 'offset after level record:',f.tell()
        iOct = iHOLL[Lev] - 1
        nLevel = iNOLL[Lev]
        nLevCells = nLevel * nchild
        ntot = ntot + nLevel

        # Skip all the oct hierarchy data
        ns = peek_record_size(f, endian='>')
        size = struct.calcsize('>i') + ns + struct.calcsize('>i')
        f.seek(f.tell()+size * nLevel)

        level_child_offsets.append(f.tell())
        # Skip the child vars data
        ns = peek_record_size(f, endian='>')
        size = struct.calcsize('>i') + ns + struct.calcsize('>i')
        f.seek(f.tell()+size * nLevel*nchild)

        # find nhydrovars
        nhydrovars = 8+2
    f.seek(offset)
    return nhydrovars, iNOLL, level_oct_offsets, level_child_offsets


def _read_art_level_info(f, level_oct_offsets, level, coarse_grid=128,
                         ncell0=None, root_level=None):
    pos = f.tell()
    f.seek(level_oct_offsets[level])
    # Get the info for this level, skip the rest
    junk, nLevel, iOct = read_vector(f, 'i', '>')

    # fortran indices start at 1

    # Skip all the oct hierarchy data
    le = np.zeros((nLevel, 3), dtype='int64')
    fl = np.ones((nLevel, 6), dtype='int64')
    iocts = np.zeros(nLevel+1, dtype='int64')
    idxa, idxb = 0, 0
    chunk = long(1e6)  # this is ~111MB for 15 dimensional 64 bit arrays
    left = nLevel
    while left > 0:
        this_chunk = min(chunk, left)
        idxb = idxa+this_chunk
        data = np.fromfile(f, dtype='>i', count=this_chunk*15)
        data = data.reshape(this_chunk, 15)
        left -= this_chunk
        le[idxa:idxb, :] = data[:, 1:4]
        fl[idxa:idxb, 1] = np.arange(idxa, idxb)
        # pad byte is last, LL2, then ioct right before it
        iocts[idxa:idxb] = data[:, -3]
        idxa = idxa+this_chunk
    del data

    # emulate fortran code
    #     do ic1 = 1 , nLevel
    #       read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
    #&                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct),
    #&                iOctLL2(iOct)
    #       iOct = iOctLL1(iOct)

    # ioct always represents the index of the next variable
    # not the current, so shift forward one index
    # the last index isn't used
    ioctso = iocts.copy()
    iocts[1:] = iocts[:-1]  # shift
    iocts = iocts[:nLevel]  # chop off the last, unused, index
    iocts[0] = iOct  # starting value

    # now correct iocts for fortran indices start @ 1
    iocts = iocts-1

    assert np.unique(iocts).shape[0] == nLevel

    # left edges are expressed as if they were on
    # level 15, so no matter what level max(le)=2**15
    # correct to the yt convention
    # le = le/2**(root_level-1-level)-1

    # try to find the root_level first
    def cfc(root_level, level, le):
        d_x = 1.0/(2.0**(root_level-level+1))
        fc = (d_x * le) - 2**(level-1)
        return fc
    if root_level is None:
        root_level = np.floor(np.log2(le.max()*1.0/coarse_grid))
        root_level = root_level.astype('int64')
        for i in range(10):
            fc = cfc(root_level, level, le)
            go = np.diff(np.unique(fc)).min() < 1.1
            if go:
                break
            root_level += 1
    else:
        fc = cfc(root_level, level, le)
    unitary_center = fc/(coarse_grid*2.0**(level-1))
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
              #[-1, -1,  1,  1, -1, -1,  1,  1],
              #[-1, -1, -1, -1,  1,  1,  1,  1]]
    # idelta = np.array(idelta)
    # if ncell0 is None:
        # ncell0 = coarse_grid**3
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


def read_particles(file, Nrow, dd=1.0, idxa=None, idxb=None):
    words = 6  # words (reals) per particle: x,y,z,vx,vy,vz
    real_size = 4  # for file_particle_data; not always true?
    np_per_page = Nrow**2  # defined in ART a_setup.h
    num_pages = os.path.getsize(file)/(real_size*words*np_per_page)

    f = np.fromfile(file, dtype='>f4').astype('float32')  # direct access
    pages = np.vsplit(np.reshape(f, (
        num_pages, words, np_per_page)), num_pages)
    data = np.squeeze(np.dstack(pages)).T  # x,y,z,vx,vy,vz
    return data[idxa:idxb, 0:3]/dd, data[idxa:idxb, 3:]


def read_star_field(file, field=None):
    data = {}
    with open(file, 'rb') as fh:
        for dtype, variables in star_struct:
            found = field in variables or field == variables
            if found:
                data[field] = read_vector(fh, dtype[1], dtype[0])
            else:
                skip(fh, endian='>')
    return data.pop(field)


def _read_child_mask_level(f, level_child_offsets, level, nLevel, nhydro_vars):
    f.seek(level_child_offsets[level])
    nvals = nLevel * (nhydro_vars + 6)  # 2 vars, 2 pads
    ioctch = np.zeros(nLevel, dtype='uint8')
    idc = np.zeros(nLevel, dtype='int32')

    chunk = long(1e6)
    left = nLevel
    width = nhydro_vars+6
    a, b = 0, 0
    while left > 0:
        chunk = min(chunk, left)
        b += chunk
        arr = np.fromfile(f, dtype='>i', count=chunk*width)
        arr = arr.reshape((width, chunk), order="F")
        assert np.all(arr[0, :] == arr[-1, :])  # pads must be equal
        idc[a:b] = arr[1, :]-1  # fix fortran indexing
        ioctch[a:b] = arr[
            2, :] == 0  # if it is above zero, then refined available
        # zero in the mask means there is refinement available
        a = b
        left -= chunk
    assert left == 0
    return idc, ioctch

nchem = 8+2
dtyp = np.dtype(">i4,>i8,>i8"+",>%sf4" % (nchem) +
                ",>%sf4" % (2)+",>i4")


def _read_child_level(
    f, level_child_offsets, level_oct_offsets, level_info, level,
    fields, domain_dimensions, ncell0, nhydro_vars=10, nchild=8,
        noct_range=None):
    # emulate the fortran code for reading cell data
    # read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar),
    #    &                 (var(i,idc), i=2,3)
    # contiguous 8-cell sections are for the same oct;
    # ie, we don't write out just the 0 cells, then the 1 cells
    # optionally, we only read noct_range to save memory
    left_index, fl, octs, nocts, root_level = _read_art_level_info(f,
                                                                   level_oct_offsets, level, coarse_grid=domain_dimensions[0])
    if noct_range is None:
        nocts = level_info[level]
        ncells = nocts*8
        f.seek(level_child_offsets[level])
        arr = np.fromfile(f, dtype=hydro_struct, count=ncells)
        assert np.all(arr['pad1'] == arr['pad2'])  # pads must be equal
        # idc = np.argsort(arr['idc']) #correct fortran indices
        # translate idc into icell, and then to iOct
        icell = (arr['idc'] >> 3) << 3
        iocts = (icell-ncell0)/nchild  # without a F correction, theres a +1
        # assert that the children are read in the same order as the octs
        assert np.all(octs == iocts[::nchild])
    else:
        start, end = noct_range
        nocts = min(end-start, level_info[level])
        end = start + nocts
        ncells = nocts*8
        skip = np.dtype(hydro_struct).itemsize*start*8
        f.seek(level_child_offsets[level]+skip)
        arr = np.fromfile(f, dtype=hydro_struct, count=ncells)
        assert np.all(arr['pad1'] == arr['pad2'])  # pads must be equal
    source = {}
    for field in fields:
        sh = (nocts, 8)
        source[field] = np.reshape(arr[field], sh, order='C').astype('float64')
    return source


def _read_root_level(f, level_offsets, level_info, nhydro_vars=10):
    nocts = level_info[0]
    f.seek(level_offsets[0])  # Ditch the header
    hvar = read_vector(f, 'f', '>')
    var = read_vector(f, 'f', '>')
    hvar = hvar.reshape((nhydro_vars, nocts*8), order="F")
    var = var.reshape((2, nocts*8), order="F")
    arr = np.concatenate((hvar, var))
    return arr

# All of these functions are to convert from hydro time var to
# proper time
sqrt = np.sqrt
sign = np.sign


def find_root(f, a, b, tol=1e-6):
    c = (a+b)/2.0
    last = -np.inf
    assert(sign(f(a)) != sign(f(b)))
    while np.abs(f(c)-last) > tol:
        last = f(c)
        if sign(last) == sign(f(b)):
            b = c
        else:
            a = c
        c = (a+b)/2.0
    return c


def quad(fintegrand, xmin, xmax, n=1e4):
    spacings = np.logspace(np.log10(xmin), np.log10(xmax), n)
    integrand_arr = fintegrand(spacings)
    val = np.trapz(integrand_arr, dx=np.diff(spacings))
    return val


def a2b(at, Om0=0.27, Oml0=0.73, h=0.700):
    def f_a2b(x):
        val = 0.5*sqrt(Om0) / x**3.0
        val /= sqrt(Om0/x**3.0 + Oml0 + (1.0 - Om0-Oml0)/x**2.0)
        return val
    # val, err = si.quad(f_a2b,1,at)
    val = quad(f_a2b, 1, at)
    return val


def b2a(bt, **kwargs):
    # converts code time into expansion factor
    # if Om0 ==1and OmL == 0 then b2a is (1 / (1-td))**2
    # if bt < -190.0 or bt > -.10:  raise 'bt outside of range'
    f_b2a = lambda at: a2b(at, **kwargs)-bt
    return find_root(f_b2a, 1e-4, 1.1)
    # return so.brenth(f_b2a,1e-4,1.1)
    # return brent.brent(f_b2a)


def a2t(at, Om0=0.27, Oml0=0.73, h=0.700):
    integrand = lambda x: 1./(x*sqrt(Oml0+Om0*x**-3.0))
    # current_time,err = si.quad(integrand,0.0,at,epsabs=1e-6,epsrel=1e-6)
    current_time = quad(integrand, 1e-4, at)
    # spacings = np.logspace(-5,np.log10(at),1e5)
    # integrand_arr = integrand(spacings)
    # current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
    current_time *= 9.779/h
    return current_time


def b2t(tb, n=1e2, logger=None, **kwargs):
    tb = np.array(tb)
    if isinstance(tb, type(1.1)):
        return a2t(b2a(tb))
    if tb.shape == ():
        return a2t(b2a(tb))
    if len(tb) < n:
        n = len(tb)
    age_min = a2t(b2a(tb.max(), **kwargs), **kwargs)
    age_max = a2t(b2a(tb.min(), **kwargs), **kwargs)
    tbs = -1.*np.logspace(np.log10(-tb.min()),
                          np.log10(-tb.max()), n)
    ages = []
    for i, tbi in enumerate(tbs):
        ages += a2t(b2a(tbi)),
        if logger:
            logger(i)
    ages = np.array(ages)
    return tbs, ages
