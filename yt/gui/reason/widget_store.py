"""
This is a place to store widgets, and to create them.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.mods import *
import weakref
from .bottle_mods import PayloadHandler, lockit
from .widget_builders import RenderingScene, get_corners, get_isocontour
from yt.visualization.plot_window import PWViewerExtJS
from yt.visualization.volume_rendering.create_spline import create_spline
import uuid

_phase_plot_mds = """<pre>
Field X: %s
Field Y: %s
Field Z: %s
Weight : %s
</pre>
"""

class WidgetStore(dict):
    def __init__(self, global_token = ""):
        self._global_token = global_token
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

    def create_slice(self, ds, center, axis, field, onmax):
        if onmax: 
            center = ds.h.find_max('Density')[1]
        else:
            center = np.array(center)
        axis = ds.coordinates.axis_id[axis.lower()]
        coord = center[axis]
        sl = ds.slice(axis, coord, center = center)
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        DLE, DRE = ds.domain_left_edge, ds.domain_right_edge
        pw = PWViewerExtJS(sl, (DLE[xax], DRE[xax], DLE[yax], DRE[yax]), 
                           setup = False, plot_type='SlicePlot')
        pw.set_current_field(field)
        field_list = list(set(ds.field_list + ds.derived_field_list))
        field_list = [dict(text = f) for f in sorted(field_list)]
        cb = pw._get_cbar_image()
        trans = pw._field_transform[pw._current_field].name
        widget_data = {'fields': field_list,
                         'initial_field': field,
                         'title': "%s Slice" % (ds),
                         'colorbar': cb,
                         'initial_transform' : trans}
        self._add_widget(pw, widget_data)

    def create_proj(self, ds, axis, field, weight):
        if weight == "None": weight = None
        axis = ds.coordinates.axis_id[axis.lower()]
        proj = ds.proj(field, axis, weight_field=weight)
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        DLE, DRE = ds.domain_left_edge, ds.domain_right_edge
        pw = PWViewerExtJS(proj, (DLE[xax], DRE[xax], DLE[yax], DRE[yax]),
                           setup = False, plot_type='ProjectionPlot')
        pw.set_current_field(field)
        field_list = list(set(ds.field_list + ds.derived_field_list))
        field_list = [dict(text = f) for f in sorted(field_list)]
        cb = pw._get_cbar_image()
        widget_data = {'fields': field_list,
                       'initial_field': field,
                       'title': "%s Projection" % (ds),
                       'colorbar': cb}
        self._add_widget(pw, widget_data)

    def create_grid_dataview(self, ds):
        levels = ds.grid_levels
        left_edge = ds.grid_left_edge
        right_edge = ds.grid_right_edge
        dimensions = ds.grid_dimensions
        cell_counts = ds.grid_dimensions.prod(axis=1)
        # This is annoying, and not ... that happy for memory.
        i = ds.index.grids[0]._id_offset
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

    def create_ds_display(self, ds):
        widget = ParameterFileWidget(ds)
        widget_data = {'fields': widget._field_list(),
                       'level_stats': widget._level_stats(),
                       'ds_info': widget._ds_info(),
                      }
        self._add_widget(widget, widget_data)

    def create_mapview(self, widget_name):
        widget = self[widget_name]
        # We want multiple maps simultaneously
        uu = "/%s/%s" % (self._global_token,
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

    def create_scene(self, ds):
        '''Creates 3D XTK-based scene'''
        widget = SceneWidget(ds)
        field_list = list(set(ds.field_list
                            + ds.derived_field_list))
        field_list.sort()
        widget_data = {'title':'Scene for %s' % ds,
                       'fields': field_list}
        self._add_widget(widget, widget_data)


class ParameterFileWidget(object):
    _ext_widget_id = None
    _widget_name = "parameterfile"

    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    def _field_list(self):
        field_list = list(set(self.ds.field_list
                            + self.ds.derived_field_list))
        field_list.sort()
        return [dict(text = field) for field in field_list]

    def _level_stats(self):
        level_data = []
        level_stats = self.ds.h.level_stats
        ngrids = float(level_stats['numgrids'].sum())
        ncells = float(level_stats['numcells'].sum())
        for level in range(self.ds.h.max_level + 1):
            cell_count = level_stats['numcells'][level]
            grid_count = level_stats['numgrids'][level]
            level_data.append({'level' : level,
                               'cell_count': int(cell_count),
                               'grid_count': int(grid_count),
                               'cell_rel': int(100*cell_count/ncells),
                               'grid_rel': int(100*grid_count/ngrids)})
        return level_data

    def _ds_info(self):
        tr = {}
        for k, v in self.ds._mrep._attrs.items():
            if isinstance(v, np.ndarray):
                tr[k] = v.tolist()
            else:
                tr[k] = v
        return tr

    def deliver_field(self, field):
        ph = PayloadHandler()
        ph.widget_payload(self,
            {'ptype':'field_info',
             'field_source': self.ds.field_info[field].get_source() })
        return

class SceneWidget(object):
    _ext_widget_id = None
    _widget_name = "scene"
    _rendering_scene = None

    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    def add_volume_rendering(self):
        return None
        self._rendering_scene = RenderingScene(self.ds, None, None)

    def deliver_rendering(self, scene_config):
        ph = PayloadHandler()
        ph.widget_payload(self, {'ptype':'rendering',
                                 'image':self._rendering_scene.snapshot()})
        return

    def deliver_isocontour(self, field, value, rel_val = False):
        ph = PayloadHandler()
        vert = get_isocontour(self.ds, field, value, rel_val)
        normals = np.empty(vert.shape)
        for i in xrange(vert.shape[0]/3):
            n = np.cross(vert[i*3,:], vert[i*3+1,:])
            normals[i*3:i*3+3,:] = n[None,:]
        ph.widget_payload(self, {'ptype':'isocontour',
                                 'binary': ['vert', 'normals'],
                                 'vert':vert,
                                 'normals':normals})

    def deliver_gridlines(self):
        ph = PayloadHandler()
        corners, levels = get_corners(self.ds)
        ph.widget_payload(self, {'ptype':'grid_lines',
                                 'binary': ['corners', 'levels'],
                                 'corners': corners,
                                 'levels': levels,
                                 'max_level': int(self.ds.h.max_level)})
        return

    def render_path(self, views, times, N):
        # Assume that path comes in as a list of matrice
        # Assume original vector is (0., 0., 1.), up is (0., 1., 0.)
        
        views = [np.array(view).transpose() for view in views]

        times = np.linspace(0.0,1.0,len(times))
                
        # This is wrong.
        reflect = np.array([[1,0,0],[0,1,0],[0,0,-1]])

        rots = np.array([R[0:3,0:3] for R in views])

        rots = np.array([np.dot(reflect,rot) for rot in rots])

        centers = np.array([np.dot(rot,R[0:3,3]) for R,rot in zip(views,rots)])

        ups = np.array([np.dot(rot,R[0:3,1]) for R,rot in zip(views,rots)])

        #print 'views'
        #for view in views: print view
        #print 'rots'
        #for rot in rots: print rot
        #print 'centers'
        #for center in centers: print center
        #print 'ups'
        #for up in ups: print up

        pos = np.empty((N,3), dtype="float64")
        uv = np.empty((N,3), dtype="float64")
        f = np.zeros((N,3), dtype="float64")
        for i in range(3):
            pos[:,i] = create_spline(times, centers[:,i], np.linspace(0.0,1.0,N))
            uv[:,i] = create_spline(times, ups[:,i], np.linspace(0.0,1.0,N))
    
        path = [pos.tolist(), f.tolist(), uv.tolist()]
    
        ph = PayloadHandler()
        ph.widget_payload(self, {'ptype':'camerapath',
                                 'data':path,})

        return


    def deliver_streamlines(self):
        ds = PayloadHandler()
        pass

