"""
ARTIO-specific IO




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler


class IOHandlerARTIO(BaseIOHandler):
    _dataset_type = "artio"

    def _read_fluid_selection(self, chunks, selector, fields):
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        cp = 0
        for onechunk in chunks:
            for artchunk in onechunk.objs:
                rv = artchunk.fill(fields, selector)
                for f in fields:
                    tr[f].resize(cp+artchunk.data_size)
                    tr[f][cp:cp+artchunk.data_size] = rv.pop(f)
                cp += artchunk.data_size
        return tr

    def _read_particle_coords(self, chunks, ptf):
        pn = "POSITION_%s"
        chunks = list(chunks)
        fields = [(ptype, pn % ax)
                  for ptype, field_list in ptf.items()
                  for ax in 'XYZ']
        for chunk in chunks: # These should be organized by grid filename
            for subset in chunk.objs:
                rv = dict(**subset.fill_particles(fields))
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (np.asarray(rv[ptype][pn % ax], dtype="=f8")
                               for ax in 'XYZ')
                    yield ptype, (x, y, z)
                    rv.pop(ptype)

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = "POSITION_%s"
        chunks = list(chunks)
        fields = [(ptype, fname) for ptype, field_list in ptf.items()
                                 for fname in field_list]
        for ptype, field_list in sorted(ptf.items()):
            for ax in 'XYZ':
                if pn % ax not in field_list:
                    fields.append((ptype, pn % ax))
        for chunk in chunks: # These should be organized by grid filename
            for subset in chunk.objs:
                rv = dict(**subset.fill_particles(fields))
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (np.asarray(rv[ptype][pn % ax], dtype="=f8")
                               for ax in 'XYZ')
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None: continue
                    for field in field_list:
                        data = np.asarray(rv[ptype][field], "=f8")
                        yield (ptype, field), data[mask]
                    rv.pop(ptype)
