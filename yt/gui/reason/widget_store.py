"""
This is a place to store widgets, and to create them.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

from yt.mods import *
import weakref
from .bottle_mods import PayloadHandler, lockit
from .widget_builders import RenderingScene 
from yt.visualization.plot_window import PWViewerExtJS
import uuid

_phase_plot_mds = """<pre>
Field X: %s
Field Y: %s
Field Z: %s
Weight : %s
</pre>
"""

class WidgetStore(dict):
    def __init__(self, repl):
        self.repl = weakref.proxy(repl)
        self.payload_handler = PayloadHandler()
        super(WidgetStore, self).__init__()

    def _add_widget(self, widget, widget_data = None):
        # We need to make sure that we aren't running in advance of a new
        # object being added.
        #uu = str(uuid.uuid1()).replace("-","_")
        #varname = "%s_%s" % (widget._widget_name, uu)
        varname = "%s_%s" % (widget._widget_name, len(self))
        widget._ext_widget_id = varname
        payload = {'type': 'widget',
                   'widget_type': widget._widget_name,
                   'varname': varname}
        widget._ext_widget_id = varname
        payload['data'] = widget_data
        self[varname] = widget
        self.payload_handler.add_payload(payload)

    def activate_multicast(self, widget_id, multicast_session, multicast_token):
        # Here we have to conduct the handshake between us and the GAE
        self.payload_handler.multicast_ids[widget_id] = (
            multicast_session, multicast_token)
        mylog.info("Multicasting %s to %s", widget_id, multicast_session)

    def create_slice(self, pf, center, axis, field, onmax):
        if onmax: 
            center = pf.h.find_max('Density')[1]
        else:
            center = na.array(center)
        axis = inv_axis_names[axis.lower()]
        coord = center[axis]
        sl = pf.h.slice(axis, coord, center = center, periodic = True)
        xax, yax = x_dict[axis], y_dict[axis]
        DLE, DRE = pf.domain_left_edge, pf.domain_right_edge
        pw = PWViewerExtJS(sl, (DLE[xax], DRE[xax], DLE[yax], DRE[yax]), setup = False)
        pw.set_current_field(field)
        field_list = list(set(pf.h.field_list + pf.h.derived_field_list))
        field_list = [dict(text = f) for f in sorted(field_list)]
        cb = pw._get_cbar_image()
        trans = pw._field_transform[pw._current_field].name
        widget_data = {'fields': field_list,
                         'initial_field': field,
                         'title': "%s Slice" % (pf),
                         'colorbar': cb,
                         'initial_transform' : trans}
        self._add_widget(pw, widget_data)

    def create_proj(self, pf, axis, field, weight):
        if weight == "None": weight = None
        axis = inv_axis_names[axis.lower()]
        proj = pf.h.proj(axis,field, weight_field=weight, periodic = True)
        xax, yax = x_dict[axis], y_dict[axis]
        DLE, DRE = pf.domain_left_edge, pf.domain_right_edge
        pw = PWViewerExtJS(proj, (DLE[xax], DRE[xax], DLE[yax], DRE[yax]),
                           setup = False)
        pw.set_current_field(field)
        field_list = list(set(pf.h.field_list + pf.h.derived_field_list))
        field_list = [dict(text = f) for f in sorted(field_list)]
        cb = pw._get_cbar_image()
        widget_data = {'fields': field_list,
                       'initial_field': field,
                       'title': "%s Projection" % (pf),
                       'colorbar': cb}
        self._add_widget(pw, widget_data)

    def create_grid_dataview(self, pf):
        levels = pf.h.grid_levels
        left_edge = pf.h.grid_left_edge
        right_edge = pf.h.grid_right_edge
        dimensions = pf.h.grid_dimensions
        cell_counts = pf.h.grid_dimensions.prod(axis=1)
        # This is annoying, and not ... that happy for memory.
        i = pf.h.grids[0]._id_offset
        vals = []
        for i, (L, LE, RE, dim, cell) in enumerate(zip(
            levels, left_edge, right_edge, dimensions, cell_counts)):
            vals.append([ int(i), int(L[0]),
                          float(LE[0]), float(LE[1]), float(LE[2]),
                          float(RE[0]), float(RE[1]), float(RE[2]),
                          int(dim[0]), int(dim[1]), int(dim[2]),
                          int(cell)] )
        varname = "gg_%s" % (len(self))
        self[varname] = None
        payload = {'type': 'widget',
                   'widget_type': 'grid_data',
                   'varname': varname, # Is just "None"
                   'data': dict(gridvals = vals),
                   }
        self.payload_handler.add_payload(payload)

    def create_pf_display(self, pf):
        widget = ParameterFileWidget(pf)
        widget_data = {'fields': widget._field_list(),
                       'level_stats': widget._level_stats(),
                       'pf_info': widget._pf_info(),
                      }
        self._add_widget(widget, widget_data)

    def create_mapview(self, widget_name):
        widget = self[widget_name]
        # We want multiple maps simultaneously
        uu = "/%s/%s" % (getattr(self.repl, "_global_token", ""),
                        str(uuid.uuid1()).replace("-","_"))
        from .pannable_map import PannableMapServer
        data = widget.data_source
        field_name = widget._current_field
        pm = PannableMapServer(data, field_name, route_prefix = uu)
        widget_data = {'prefix': uu, 'field':field_name}
        self._add_widget(pm, widget_data)

    def create_phase(self, obj, field_x, field_y, field_z, weight):
        from yt.visualization.profile_plotter import PhasePlotterExtWidget
        pp = PhasePlotterExtWidget(obj, field_x, field_y, field_z, weight)
        mds = _phase_plot_mds % (field_x, field_y, field_z, 
                                 pp._initial_weight)
        widget_data = {'title': "%s Phase Plot" % (obj),
                       'metadata_string': mds}
        self._add_widget(pp, widget_data)

    def create_scene(self, pf):
        '''Creates 3D XTK-based scene'''
        widget = SceneWidget(pf)
        widget_data = {'title':'Scene for %s' % pf}
        self._add_widget(widget, widget_data)


class ParameterFileWidget(object):
    _ext_widget_id = None
    _widget_name = "parameterfile"

    def __init__(self, pf):
        self.pf = weakref.proxy(pf)

    def _field_list(self):
        field_list = list(set(self.pf.h.field_list
                            + self.pf.h.derived_field_list))
        field_list.sort()
        return [dict(text = field) for field in field_list]

    def _level_stats(self):
        level_data = []
        level_stats = self.pf.h.level_stats
        ngrids = float(level_stats['numgrids'].sum())
        ncells = float(level_stats['numcells'].sum())
        for level in range(self.pf.h.max_level + 1):
            cell_count = level_stats['numcells'][level]
            grid_count = level_stats['numgrids'][level]
            level_data.append({'level' : level,
                               'cell_count': int(cell_count),
                               'grid_count': int(grid_count),
                               'cell_rel': int(100*cell_count/ncells),
                               'grid_rel': int(100*grid_count/ngrids)})
        return level_data

    def _pf_info(self):
        tr = {}
        for k, v in self.pf._mrep._attrs.items():
            if isinstance(v, na.ndarray):
                tr[k] = v.tolist()
            else:
                tr[k] = v
        return tr

    def deliver_field(self, field):
        ph = PayloadHandler()
        ph.widget_payload(self,
            {'ptype':'field_info',
             'field_source': self.pf.field_info[field].get_source() })
        return

class SceneWidget(object):
    _ext_widget_id = None
    _widget_name = "scene"
    _rendering_scene = None

    def __init__(self, pf):
        self.pf = weakref.proxy(pf)

    def add_volume_rendering(self):
        return None
        self._rendering_scene = RenderingScene(self.pf, None, None)

    def deliver_rendering(self, scene_config):
        ph = PayloadHandler()
        ph.widget_payload(self, {'ptype':'rendering',
                                 'image':self._rendering_scene.snapshot()})
        return

    def deliver_gridlines(self):
        ph = PayloadHandler()
        ph.widget_payload(self, {'ptype':'grid_lines',
                                 'vertices': get_grid_vertices(),
                                 'colors': get_grid_colors()})
        return

    def deliver_streamlines(self):
        pf = PayloadHandler()
        pass

