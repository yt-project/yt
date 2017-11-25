"""
RAMSES-specific IO



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import defaultdict
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.physical_ratios import cm_per_km, cm_per_mpc
import yt.utilities.fortran_utils as fpu
from yt.utilities.lib.cosmology_time import \
    get_ramses_ages
from yt.utilities.exceptions import YTFieldTypeNotFound, YTParticleOutputFormatNotImplemented, \
    YTFileNotParseable
from yt.extern.six import PY3
import re

if PY3:
    from io import BytesIO as IO
else:
    from cStringIO import StringIO as IO

def _ramses_particle_file_handler(fname, foffsets, data_types,
                                  subset, fields):
    '''General file handler, called by _read_particle_subset

    Parameters
    ----------
    fname : string
        filename to read from
    foffsets: dict
        Offsets in file of the fields
    data_types: dict
         Data type of the fields
    subset: ``RAMSESDomainSubset``
         A RAMSES domain subset object
    fields: list of tuple
         The fields to read
    '''
    tr = {}
    with open(fname, "rb") as f:
        # We do *all* conversion into boxlen here.
        # This means that no other conversions need to be applied to convert
        # positions into the same domain as the octs themselves.
        for field in sorted(fields, key=lambda a: foffsets[a]):
            f.seek(foffsets[field])
            dt = data_types[field]
            tr[field] = fpu.read_vector(f, dt)
            if field[1].startswith("particle_position"):
                np.divide(tr[field], subset.domain.ds["boxlen"], tr[field])
            cosmo = subset.domain.ds.cosmological_simulation
            if cosmo == 1 and field[1] == "particle_age":
                tf = subset.domain.ds.t_frw
                dtau = subset.domain.ds.dtau
                tauf = subset.domain.ds.tau_frw
                tsim = subset.domain.ds.time_simu
                h100 = subset.domain.ds.hubble_constant
                nOver2 = subset.domain.ds.n_frw/2
                t_scale = 1./(h100 * 100 * cm_per_km / cm_per_mpc)/subset.domain.ds['unit_t']
                ages = tr[field]
                tr[field] = get_ramses_ages(tf,tauf,dtau,tsim,t_scale,ages,nOver2,len(ages))
    return tr


class IOHandlerRAMSES(BaseIOHandler):
    _dataset_type = "ramses"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Chunks in this case will have affiliated domain subset objects
        # Each domain subset will contain a hydro_offset array, which gives
        # pointers to level-by-level hydro information
        tr = defaultdict(list)
        for chunk in chunks:
            for subset in chunk.objs:
                # Now we read the entire thing
                f = open(subset.domain.hydro_fn, "rb")
                # This contains the boundary information, so we skim through
                # and pick off the right vectors
                content = IO(f.read())
                rv = subset.fill(content, fields, selector)
                for ft, f in fields:
                    d = rv.pop(f)
                    mylog.debug("Filling %s with %s (%0.3e %0.3e) (%s zones)",
                        f, d.size, d.min(), d.max(), d.size)
                    tr[(ft, f)].append(d)
        d = {}
        for field in fields:
            d[field] = np.concatenate(tr.pop(field))
        return d

    def _read_particle_coords(self, chunks, ptf):
        pn = "particle_position_%s"
        fields = [(ptype, "particle_position_%s" % ax)
                  for ptype, field_list in ptf.items()
                  for ax in 'xyz']
        for chunk in chunks:
            for subset in chunk.objs:
                rv = self._read_particle_subset(subset, fields)
                for ptype in sorted(ptf):
                    yield ptype, (rv[ptype, pn % 'x'],
                                  rv[ptype, pn % 'y'],
                                  rv[ptype, pn % 'z'])

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = "particle_position_%s"
        chunks = list(chunks)
        fields = [(ptype, fname) for ptype, field_list in ptf.items()
                                 for fname in field_list]
        for ptype, field_list in sorted(ptf.items()):
            for ax in 'xyz':
                if pn % ax not in field_list:
                    fields.append((ptype, pn % ax))
        for chunk in chunks:
            for subset in chunk.objs:
                rv = self._read_particle_subset(subset, fields)
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (np.asarray(rv[ptype, pn % ax], "=f8")
                               for ax in 'xyz')
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                       mask = []
                    for field in field_list:
                        data = np.asarray(rv.pop((ptype, field))[mask], "=f8")
                        yield (ptype, field), data


    def _read_particle_subset(self, subset, fields):
        '''Read the particle files.'''
        tr = {}

        # Sequential read depending on particle type (io or sink)
        for ptype in set(f[0] for f in fields):

            # Select relevant fiels
            subs_fields = filter(lambda f: f[0] == ptype, fields)

            if ptype == 'io':
                fname = subset.domain.part_fn
                foffsets = subset.domain.particle_field_offsets
                data_types = subset.domain.particle_field_types

            elif ptype == 'sink':
                fname = subset.domain.sink_fn
                foffsets = subset.domain.sink_field_offsets
                data_types = subset.domain.sink_field_types

            else:
                # Raise here an exception
                raise YTFieldTypeNotFound(ptype)

            tr.update(_ramses_particle_file_handler(
                fname, foffsets, data_types, subset, subs_fields))

        return tr

VERSION_RE = re.compile('# version: *(\d+)')
VAR_DESC_RE = re.compile(r'\s*(\d+),\s*(\w+),\s*(\w+)')
def _read_part_file_descriptor(fname):
    """
    Read the particle file descriptor and returns the array of the fields found.
    """
    # Mapping
    mapping = [
        ('position_x', 'particle_position_x'),
        ('position_y', 'particle_position_y'),
        ('position_z', 'particle_position_z'),
        ('velocity_x', 'particle_velocity_x'),
        ('velocity_y', 'particle_velocity_y'),
        ('velocity_z', 'particle_velocity_z'),
        ('mass', 'particle_mass'),
        ('identity', 'particle_identity'),
        ('levelp', 'particle_level'),
        ('family', 'particle_family'),
        ('tag', 'particle_tag')
    ]
    # Convert in dictionary
    mapping = {k: v for k, v in mapping}

    with open(fname, 'r') as f:
        f = open(fname, 'r')
        line = f.readline()
        tmp = VERSION_RE.match(line)
        mylog.info('Reading part file descriptor.')
        if not tmp:
            raise YTParticleOutputFormatNotImplemented()

        version = int(tmp.group(1))

        if version == 1:
            # Skip one line (containing the headers)
            line = f.readline()
            fields = []
            for i, line in enumerate(f.readlines()):
                tmp = VAR_DESC_RE.match(line)
                if not tmp:
                    raise YTFileNotParseable(fname, i+1)

                # ivar = tmp.group(1)
                varname = tmp.group(2)
                dtype = tmp.group(3)

                if varname in mapping:
                    varname = mapping[varname]
                else:
                    varname = 'particle_%s' % varname

                fields.append((varname, dtype))
        else:
            raise YTParticleOutputFormatNotImplemented()

    return fields
