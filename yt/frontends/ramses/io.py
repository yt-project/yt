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
from yt.extern.six import PY3

if PY3:
    from io import BytesIO as IO
else:
    from cStringIO import StringIO as IO

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
        f = open(subset.domain.part_fn, "rb")
        foffsets = subset.domain.particle_field_offsets
        tr = {}
        # We do *all* conversion into boxlen here.
        # This means that no other conversions need to be applied to convert
        # positions into the same domain as the octs themselves.
        for field in sorted(fields, key = lambda a: foffsets[a]):
            f.seek(foffsets[field])
            dt = subset.domain.particle_field_types[field]
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
