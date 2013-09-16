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
import yt.utilities.fortran_utils as fpu
import cStringIO

class IOHandlerRAMSES(BaseIOHandler):
    _data_style = "ramses"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # Chunks in this case will have affiliated domain subset objects
        # Each domain subset will contain a hydro_offset array, which gives
        # pointers to level-by-level hydro information
        tr = defaultdict(list)
        cp = 0
        for chunk in chunks:
            for subset in chunk.objs:
                # Now we read the entire thing
                f = open(subset.domain.hydro_fn, "rb")
                # This contains the boundary information, so we skim through
                # and pick off the right vectors
                content = cStringIO.StringIO(f.read())
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

    def _read_particle_selection(self, chunks, selector, fields):
        size = 0
        masks = {}
        chunks = list(chunks)
        pos_fields = [("all","particle_position_%s" % ax) for ax in "xyz"]
        for chunk in chunks:
            for subset in chunk.objs:
                # We read the whole thing, then feed it back to the selector
                selection = self._read_particle_subset(subset, pos_fields)
                mask = selector.select_points(
                    selection["all", "particle_position_x"],
                    selection["all", "particle_position_y"],
                    selection["all", "particle_position_z"])
                if mask is None: continue
                #print "MASK", mask
                size += mask.sum()
                masks[id(subset)] = mask
        # Now our second pass
        tr = {}
        pos = 0
        for chunk in chunks:
            for subset in chunk.objs:
                selection = self._read_particle_subset(subset, fields)
                mask = masks.pop(id(subset), None)
                if mask is None: continue
                count = mask.sum()
                for field in fields:
                    ti = selection.pop(field)[mask]
                    if field not in tr:
                        dt = subset.domain.particle_field_types[field]
                        tr[field] = np.empty(size, dt)
                    tr[field][pos:pos+count] = ti
                pos += count
        return tr

    def _read_particle_subset(self, subset, fields):
        f = open(subset.domain.part_fn, "rb")
        foffsets = subset.domain.particle_field_offsets
        tr = {}
        # We do *all* conversion into boxlen here.
        # This means that no other conversions need to be applied to convert
        # positions into the same domain as the octs themselves.
        for field in fields:
            f.seek(foffsets[field])
            dt = subset.domain.particle_field_types[field]
            tr[field] = fpu.read_vector(f, dt)
            if field[1].startswith("particle_position"):
                np.divide(tr[field], subset.domain.pf["boxlen"], tr[field])
        return tr
