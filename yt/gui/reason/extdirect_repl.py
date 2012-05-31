"""
A read-eval-print-loop that is served up through Bottle and accepts its
commands through ExtDirect calls

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

import json
import os
import cStringIO
import logging
import uuid
import numpy as na
import time
import urllib
import urllib2
import pprint
import traceback
import tempfile
import base64
import imp
import threading
import Queue
import zipfile

from yt.funcs import *
from yt.utilities.logger import ytLogger, ufstring
from yt.utilities.definitions import inv_axis_names
from yt.visualization.image_writer import apply_colormap
from yt.visualization.api import Streamlines

from .bottle_mods import preroute, BottleDirectRouter, notify_route, \
                         PayloadHandler
from yt.utilities.bottle import response, request, route, static_file
from .basic_repl import ProgrammaticREPL

try:
    import pygments
    import pygments.lexers
    import pygments.formatters
    def _highlighter():
        pl = pygments.lexers.PythonLexer()
        hf = pygments.formatters.HtmlFormatter(
                linenos='table', linenospecial=2)
        def __highlighter(a):
            return pygments.highlight(a, pl, hf)
        return __highlighter, hf.get_style_defs()
    # We could add an additional '.highlight_pyg' in the call
    highlighter, highlighter_css = _highlighter()
except ImportError:
    highlighter = lambda a: a
    highlight_css = ''

local_dir = os.path.dirname(__file__)

class MethodLock(object):
    _shared_state = {}
    locks = None

    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        if self.locks is None: self.locks = {}

    def __call__(self, func):
        if str(func) not in self.locks:
            self.locks[str(func)] = threading.Lock()
        @wraps(func)
        def locker(*args, **kwargs):
            print "Acquiring lock on %s" % (str(func))
            with self.locks[str(func)]:
                rv = func(*args, **kwargs)
            print "Regained lock on %s" % (str(func))
            return rv
        return locker

lockit = MethodLock()

class ExecutionThread(threading.Thread):
    def __init__(self, repl):
        self.repl = repl
        self.queue = Queue.Queue()
        threading.Thread.__init__(self)
        self.daemon = True

    def run(self):
        while 1:
            #print "Checking for a queue ..."
            try:
                task = self.queue.get(True, 1)
            except Queue.Empty:
                if self.repl.stopped: return
                continue
            #print "Received the task", task
            if task['type'] == 'code':
                self.execute_one(task['code'], task['hide'])
                self.queue.task_done()
            elif task['type'] == 'add_widget':
                #print "Adding new widget"
                self.queue.task_done()
                new_code = self.repl._add_widget(
                    task['name'], task['widget_data_name'])
                #print "Got this command:", new_code
                self.execute_one(new_code, hide=True)
                #print "Executed!"

    def wait(self):
        self.queue.join()

    def execute_one(self, code, hide):
        self.repl.executed_cell_texts.append(code)
        result = ProgrammaticREPL.execute(self.repl, code)
        if self.repl.debug:
            print "==================== Cell Execution ===================="
            print code
            print "====================                ===================="
            print result
            print "========================================================"
        if not hide:
            self.repl.payload_handler.add_payload(
                {'type': 'cell',
                 'output': result,
                 'input': highlighter(code),
                 'raw_input': code},
                )

def deliver_image(im):
    if hasattr(im, 'read'):
        img_data = base64.b64encode(im.read())
    elif isinstance(im, types.StringTypes) and \
         im.endswith(".png"):
        img_data = base64.b64encode(open(im).read())
    elif isinstance(im, types.StringTypes):
        img_data = im
    else:
        raise RuntimeError
    ph = PayloadHandler()
    payload = {'type':'png_string',
               'image_data':img_data}
    ph.add_payload(payload)

def reason_pylab():
    def _canvas_deliver(canvas):
        tf = tempfile.TemporaryFile()
        canvas.print_png(tf)
        tf.seek(0)
        img_data = base64.b64encode(tf.read())
        tf.close()
        deliver_image(img_data)
    def reason_draw_if_interactive():
        if matplotlib.is_interactive():
            figManager =  Gcf.get_active()
            if figManager is not None:
                _canvas_deliver(figManager.canvas)
    def reason_show(mainloop = True):
        # We ignore mainloop here
        for manager in Gcf.get_all_fig_managers():
            _canvas_deliver(manager.canvas)
    # Matplotlib has very nice backend overriding.
    # We should really use that.  This is just a hack.
    import matplotlib
    new_agg = imp.new_module("reason_agg")
    import matplotlib.backends.backend_agg as bagg
    new_agg.__dict__.update(bagg.__dict__)
    new_agg.__dict__.update(
        {'show': reason_show,
         'draw_if_interactive': reason_draw_if_interactive})
    sys.modules["reason_agg"] = new_agg
    bagg.draw_if_interactive = reason_draw_if_interactive
    from matplotlib._pylab_helpers import Gcf
    import pylab, matplotlib
    matplotlib.rcParams["backend"] = "module://reason_agg"
    pylab.switch_backend("module://reason_agg")

class ExtDirectREPL(ProgrammaticREPL, BottleDirectRouter):
    _skip_expose = ('index')
    my_name = "ExtDirectREPL"
    timeout = 660 # a minute longer than the rocket server timeout
    server = None
    stopped = False
    debug = False
    _heartbeat_timer = None

    def __init__(self, base_extjs_path, locals=None):
        # First we do the standard initialization
        self.extjs_file = zipfile.ZipFile(os.path.join(
            base_extjs_path, "ext-4.1.0-gpl.zip"), 'r')
        ProgrammaticREPL.__init__(self, locals)
        # Now, since we want to only preroute functions we know about, and
        # since they have different arguments, and most of all because we only
        # want to add them to the routing tables (which are a singleton for the
        # entire interpreter state) we apply all the pre-routing now, rather
        # than through metaclasses or other fancy decorating.
        preroute_table = dict(index = ("/", "GET"),
                              _help_html = ("/help.html", "GET"),
                              _myapi = ("/ext-repl-api.js", "GET"),
                              _session_py = ("/session.py", "GET"),
                              _highlighter_css = ("/highlighter.css", "GET"),
                              _extjs = ("/resources/extjs-4.1.0/:path#.+#", "GET"),
                              _app = ("/:path#.+#", "GET"),
                              )
        for v, args in preroute_table.items():
            preroute(args[0], method=args[1])(getattr(self, v))
        # This has to be routed to the root directory
        self.api_url = "repl"
        BottleDirectRouter.__init__(self)
        self.pflist = ExtDirectParameterFileList()
        self.executed_cell_texts = []
        self.payload_handler = PayloadHandler()
        self.execution_thread = ExecutionThread(self)
        # Now we load up all the yt.mods stuff, but only after we've finished
        # setting up.
        reason_pylab()
        self.execute("from yt.mods import *\nimport pylab\npylab.ion()")
        self.execute("from yt.data_objects.static_output import _cached_pfs", hide = True)
        self.execute("data_objects = []", hide = True)
        self.locals['load_script'] = ext_load_script
        self.locals['deliver_image'] = deliver_image

    def activate(self):
        self._setup_logging_handlers()
        # Setup our heartbeat
        self.last_heartbeat = time.time()
        self._check_heartbeat()
        self.execution_thread.start()

    def exception_handler(self, exc):
        result = {'type': 'cell',
                  'input': 'ERROR HANDLING IN REASON',
                  'output': traceback.format_exc()}
        return result

    def _setup_logging_handlers(self):
        handler = PayloadLoggingHandler()
        formatter = logging.Formatter(ufstring)
        handler.setFormatter(formatter)
        ytLogger.addHandler(handler)

    def index(self):
        root = os.path.join(local_dir, "html")
        return static_file("index.html", root)

    def heartbeat(self):
        self.last_heartbeat = time.time()
        if self.debug: print "### Heartbeat ... started: %s" % (time.ctime())
        for i in range(30):
            # Check for stop
            if self.debug: print "    ###"
            if self.stopped: return {'type':'shutdown'} # No race condition
            if self.payload_handler.event.wait(1): # One second timeout
                if self.debug: print "    ### Delivering payloads"
                rv = self.payload_handler.deliver_payloads()
                if self.debug: print "    ### Got back, returning"
                return rv
        if self.debug: print "### Heartbeat ... finished: %s" % (time.ctime())
        return []

    def _check_heartbeat(self):
        if self.server is not None:
            if not all((s._monitor.is_alive() for s in self.server.values())):
                self.shutdown()
                return
        if time.time() - self.last_heartbeat > self.timeout:
            print "Shutting down after a timeout of %s" % (self.timeout)
            #sys.exit(0)
            # Still can't shut down yet, because bottle doesn't return the
            # server instance by default.
            self.shutdown()
            return
        if self._heartbeat_timer is not None: return
        self._heartbeat_timer = threading.Timer(10, self._check_heartbeat)
        self._heartbeat_timer.start()

    def shutdown(self):
        if self.server is None:
            return
        self._heartbeat_timer.cancel()
        self.stopped = True
        self.payload_handler.event.set()
        for v in self.server.values():
            v.stop()
        for t in threading.enumerate():
            print "Found a living thread:", t

    def _help_html(self):
        root = os.path.join(local_dir, "html")
        return static_file("help.html", root)

    def _extjs(self, path):
        pp = os.path.join("extjs-4.1.0", path)
        try:
            f = self.extjs_file.open(pp)
        except KeyError:
            response.status = 404
            return
        if path[-4:].lower() in (".png", ".gif", ".jpg"):
            response.headers['Content-Type'] = "image/%s" % (path[-3:].lower())
        elif path[-4:].lower() == ".css":
            response.headers['Content-Type'] = "text/css"
        elif path[-3:].lower() == ".js":
            response.headers['Content-Type'] = "text/javascript"
        return f.read()

    def _app(self, path):
        root = os.path.join(local_dir, "html")
        return static_file(path, root)

    def _highlighter_css(self):
        response.headers['Content-Type'] = "text/css"
        return highlighter_css

    def execute(self, code, hide = False):
            task = {'type': 'code',
                    'code': code,
                    'hide': hide}
            self.execution_thread.queue.put(task)
            return dict(status = True)

    def get_history(self):
        return self.executed_cell_texts[:]

    @lockit
    def save_session(self, filename):
        if filename.startswith('~'):
            filename = os.path.expanduser(filename)
        elif not filename.startswith('/'):
            filename = os.path.join(os.getcwd(), filename)
        if os.path.exists(filename):
            return {'status': 'FAIL', 'filename': filename,
                    'error': 'File exists!'}
        try:
            f = open(filename, 'w')
            f.write("\n######\n".join(self.executed_cell_texts))
            f.close()
        except IOError as (errno, strerror):
            return {'status': 'FAIL', 'filename': filename,
                    'error': strerror}
        except:
            return {'status': 'FAIL', 'filename': filename,
                    'error': 'Unexpected error.'}
        return {'status': 'SUCCESS', 'filename': filename}

    @lockit
    def paste_session(self):
        import xmlrpclib, cStringIO
        p = xmlrpclib.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs = cs.getvalue()
        ret = p.pastes.newPaste('python', cs, None, '', '', True)
        site = "http://paste.yt-project.org/show/%s" % ret
        return {'status': 'SUCCESS', 'site': site}

    @lockit
    def paste_text(self, to_paste):
        import xmlrpclib, cStringIO
        p = xmlrpclib.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
        ret = p.pastes.newPaste('python', to_paste, None, '', '', True)
        site = "http://paste.yt-project.org/show/%s" % ret
        return {'status': 'SUCCESS', 'site': site}

    _api_key = 'f62d550859558f28c4c214136bc797c7'
    def upload_image(self, image_data, caption):
        if not image_data.startswith("data:"): return {'uploaded':False}
        prefix = "data:image/png;base64,"
        image_data = image_data[len(prefix):]
        parameters = {'key':self._api_key, 'image':image_data, type:'base64',
                      'caption': caption, 'title': "Uploaded Image from reason"}
        data = urllib.urlencode(parameters)
        req = urllib2.Request('http://api.imgur.com/2/upload.json', data)
        try:
            response = urllib2.urlopen(req).read()
        except urllib2.HTTPError as e:
            print "ERROR", e
            return {'uploaded':False}
        rv = json.loads(response)
        rv['uploaded'] = True
        pprint.pprint(rv)
        return rv

    @lockit
    def _session_py(self):
        cs = cStringIO.StringIO()
        cs.write("\n######\n".join(self.executed_cell_texts))
        cs.seek(0)
        response.headers["content-disposition"] = "attachment;"
        return cs

    def _add_widget(self, widget_name, widget_data_name = None):
        # We need to make sure that we aren't running in advance of a new
        # object being added.
        self.execution_thread.queue.join()
        widget = self.locals[widget_name]
        uu = str(uuid.uuid1()).replace("-","_")
        varname = "%s_%s" % (widget._widget_name, uu)
        widget._ext_widget_id = varname
        # THIS BREAKS THE SCRIPT DOWNLOAD!
        # We need to make the variable be bound via an execution mechanism
        command = "%s = %s\n" % (varname, widget_name)
        payload = {'type': 'widget',
                   'widget_type': widget._widget_name,
                   'varname': varname,
                   'data': None}
        widget._ext_widget_id = varname
        if widget_data_name is not None:
            payload['data'] = self.locals[widget_data_name]
        self.payload_handler.add_payload(payload)
        return command

    @lockit
    def create_phase(self, objname, field_x, field_y, field_z, weight):
        if weight == "None": weight = None
        else: weight = "'%s'" % (weight)
        funccall = """
        _tfield_x = "%(field_x)s"
        _tfield_y = "%(field_y)s"
        _tfield_z = "%(field_z)s"
        _tweight = %(weight)s
        _tobj = %(objname)s
        _tpf = _tobj.pf
        from yt.visualization.profile_plotter import PhasePlotterExtWidget
        _tpp = PhasePlotterExtWidget(_tobj, _tfield_x, _tfield_y, _tfield_z, _tweight)
        _tfield_list = list(set(_tpf.h.field_list + _tpf.h.derived_field_list))
        _tfield_list.sort()
        _twidget_data = {'title': "%%s Phase Plot" %% (_tobj)}
        """ % dict(objname = objname, field_x = field_x, field_y = field_y,
                   field_z = field_z, weight = weight)
        funccall = "\n".join(line.strip() for line in funccall.splitlines())
        self.execute(funccall, hide=True)
        self.execution_thread.queue.put({'type': 'add_widget',
                                         'name': '_tpp',
                                         'widget_data_name': '_twidget_data'})

    @lockit
    def create_proj(self, pfname, axis, field, weight, onmax):
        if weight == "None": weight = None
        else: weight = "'%s'" % (weight)
        if not onmax:
            center_string = None
        else:
            center_string = "_tpf.h.find_max('Density')[1]"
        funccall = """
        _tpf = %(pfname)s
        _taxis = %(axis)s
        _tfield = "%(field)s"
        _tweight = %(weight)s
        _tcen = %(center_string)s
        _tsl = _tpf.h.proj(_taxis,_tfield, weight_field=_tweight, periodic = True, center=_tcen)
        _txax, _tyax = x_dict[_taxis], y_dict[_taxis]
        DLE, DRE = _tpf.domain_left_edge, _tpf.domain_right_edge
        from yt.visualization.plot_window import PWViewerExtJS
        _tpw = PWViewerExtJS(_tsl, (DLE[_txax], DRE[_txax], DLE[_tyax], DRE[_tyax]), setup = False)
        _tpw.set_current_field("%(field)s")
        _tfield_list = list(set(_tpf.h.field_list + _tpf.h.derived_field_list))
        _tfield_list.sort()
        _tcb = _tpw._get_cbar_image()
        _twidget_data = {'fields': _tfield_list,
                         'initial_field': _tfield,
                         'title': "%%s Projection" %% (_tpf),
                         'colorbar': _tcb}
        """ % dict(pfname = pfname,
                   axis = inv_axis_names[axis],
                   weight = weight,
                   center_string = center_string,
                   field=field)
        # There is a call to do this, but I have forgotten it ...
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        self.execution_thread.queue.put({'type': 'add_widget',
                                         'name': '_tpw',
                                         'widget_data_name': '_twidget_data'})

    @lockit
    def create_mapview(self, widget_name):
        # We want multiple maps simultaneously
        uu = "/%s/%s" % (getattr(self, "_global_token", ""),
                        str(uuid.uuid1()).replace("-","_"))
        from .pannable_map import PannableMapServer
        data = self.locals[widget_name].data_source
        field_name = self.locals[widget_name]._current_field
        pm = PannableMapServer(data, field_name, route_prefix = uu)
        self.locals['_tpm'] = pm
        self.locals['_twidget_data'] = {'prefix': uu, 'field':field_name}
        self.execution_thread.queue.put({'type': 'add_widget',
                                         'name': '_tpm',
                                         'widget_data_name': '_twidget_data'})

    @lockit
    def create_slice(self, pfname, center, axis, field, onmax):
        if not onmax: 
            center_string = \
              "na.array([%(c1)0.20f,%(c2)0.20f, %(c3)0.20f],dtype='float64')" \
                % dict(c1 = float(center[0]),
                       c2 = float(center[1]),
                       c3 = float(center[2]))
        else:
            center_string = "_tpf.h.find_max('Density')[1]"
        funccall = """
        _tpf = %(pfname)s
        _taxis = %(axis)s
        _tfield = "%(field)s"
        _tcenter = %(center_string)s
        _tcoord = _tcenter[_taxis]
        _tsl = _tpf.h.slice(_taxis, _tcoord, center = _tcenter, periodic = True)
        _txax, _tyax = x_dict[_taxis], y_dict[_taxis]
        DLE, DRE = _tpf.domain_left_edge, _tpf.domain_right_edge
        from yt.visualization.plot_window import PWViewerExtJS
        _tpw = PWViewerExtJS(_tsl, (DLE[_txax], DRE[_txax], DLE[_tyax], DRE[_tyax]), setup = False)
        _tpw.set_current_field("%(field)s")
        _tfield_list = list(set(_tpf.h.field_list + _tpf.h.derived_field_list))
        _tfield_list.sort()
        _tcb = _tpw._get_cbar_image()
        _ttrans = _tpw._field_transform[_tpw._current_field].name
        _twidget_data = {'fields': _tfield_list,
                         'initial_field': _tfield,
                         'title': "%%s Slice" %% (_tpf),
                         'colorbar': _tcb,
                         'initial_transform' : _ttrans}
        """ % dict(pfname = pfname,
                   center_string = center_string,
                   axis = inv_axis_names[axis],
                   field=field)
        # There is a call to do this, but I have forgotten it ...
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        self.execution_thread.queue.put({'type': 'add_widget',
                                         'name': '_tpw',
                                         'widget_data_name': '_twidget_data'})

    @lockit
    def create_isocontours(self, pfname, field, value, sampling_field):
        funccall = """
        _tpf = %(pfname)s
        _tfield = "%(field)s"
        _tvalue = %(value)s
        _tsample_values = "%(sampling_field)s"
        _tdd = _tpf.h.all_data()
        _tiso = _tdd.extract_isocontours(_tfield, _tvalue, rescale = True,
                                         sample_values = _tsample_values)
        from yt.funcs import YTEmptyClass
        _tpw = YTEmptyClass()
        print "GOT TPW"
        _tpw._widget_name = 'isocontour_viewer'
        _tpw._ext_widget_id = None
        _tverts = _tiso[0].ravel().tolist()
        _tc = (apply_colormap(na.log10(_tiso[1]))).squeeze()
        _tcolors = na.empty((_tc.shape[0] * 3, 4), dtype='float32')
        _tcolors[0::3,:] = _tc
        _tcolors[1::3,:] = _tc
        _tcolors[2::3,:] = _tc
        _tcolors = (_tcolors.ravel()/255.0).tolist()
        _twidget_data = {'vertex_positions': _tverts, 'vertex_colors': _tcolors}
        """ % dict(pfname=pfname, value=value, sampling_field=sampling_field, field=field)
        # There is a call to do this, but I have forgotten it ...
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        self.execution_thread.queue.put({'type': 'add_widget',
                                         'name' : '_tpw',
                                         'widget_data_name': '_twidget_data'})


    @lockit
    def create_grid_dataview(self, pfname):
        funccall = """
        _tpf = %(pfname)s
        """ % dict(pfname = pfname)
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        self.execution_thread.queue.join()
        pf = self.locals['_tpf']
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
        uu = str(uuid.uuid1()).replace("-","_")
        varname = "gg_%s" % (uu)
        payload = {'type': 'widget',
                   'widget_type': 'grid_data',
                   'varname': varname, # Is just "None"
                   'data': dict(gridvals = vals),
                   }
        self.execute("%s = None\n" % (varname), hide=True)
        self.payload_handler.add_payload(payload)

    @lockit
    def create_grid_viewer(self, pfname):
        funccall = """
        _tpf = %(pfname)s
        """ % dict(pfname = pfname)
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        self.execution_thread.queue.join()
        pf = self.locals['_tpf']
        corners = pf.h.grid_corners
        levels = pf.h.grid_levels
        colors = apply_colormap(levels*1.0,
                                color_bounds=[0,pf.h.max_level],
                                cmap_name="algae").repeat(24,axis=0)[:,0,:]*1.0/255.
        colors[:,3]=0.7
        colors = colors.ravel().tolist()
        
        vertices = []
        
        trans  = [0, 1, 2, 7, 5, 6, 3, 4]
        order  = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]

        for g in xrange(corners.shape[2]):
            for c in order:
                ci = trans[c]
                vertices.append(corners[ci,:,g])
        vertices = na.concatenate(vertices).tolist()
        uu = str(uuid.uuid1()).replace("-","_")
        varname = "gv_%s" % (uu)
        payload = {'type': 'widget',
                   'widget_type': 'grid_viewer',
                   'varname': varname, # Is just "None"
                   'data': dict(n_vertices = len(vertices)/3,
                                vertex_positions = vertices,
                                vertex_colors = colors)
                   }
        self.execute("%s = None\n" % (varname), hide=True)
        self.payload_handler.add_payload(payload)

    @lockit
    def create_streamline_viewer(self, pfname):
        funccall = """
        _tpf = %(pfname)s
        """ % dict(pfname = pfname)
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = True)
        pf = self.locals['_tpf']

        c = na.array([0.5]*3)
        N = 100
        scale = 1.0
        pos_dx = na.random.random((N,3))*scale-scale/2.
        pos = c+pos_dx
        
        SL = Streamlines(pf,pos,'x-velocity', 'y-velocity', 'z-velocity', length=1.0, get_magnitude=True)
        SL.integrate_through_volume()
        streamlist=[]
        stream_lengths = []
        for i,stream in enumerate(SL.streamlines):
            stream_lengths.append( stream[na.all(stream != 0.0, axis=1)].shape[0])
        streamlist = SL.streamlines.flatten()
        streamlist = streamlist[streamlist!=0.0].tolist()

        stream_colors = SL.magnitudes.flatten()
        stream_colors = na.log10(stream_colors[stream_colors > 0.0])
        stream_colors = apply_colormap(stream_colors, cmap_name='algae')
        stream_colors = stream_colors*1./255.
        stream_colors[:,:,3] = 0.8
        stream_colors = stream_colors.flatten().tolist()

        uu = str(uuid.uuid1()).replace("-","_")
        varname = "sl_%s" % (uu)
        payload = {'type': 'widget',
                   'widget_type': 'streamline_viewer',
                   'varname': varname, # Is just "None"
                   'data': dict(n_streamlines = SL.streamlines.shape[0],
                                stream_positions = streamlist,
                                stream_colors = stream_colors,
                                stream_lengths = stream_lengths)
                   }
        self.execute("%s = None\n" % (varname), hide=True)
        self.payload_handler.add_payload(payload)

    @lockit
    def object_creator(self, pfname, objtype, objargs):
        funccall = "_tobjargs = {}\n"
        for argname, argval in objargs.items():
            # These arguments may need further sanitization
            if isinstance(argval, types.StringTypes):
                argval = "'%s'" % argval
            funccall += "_tobjargs['%(argname)s'] = %(argval)s\n" % dict(
                    argname = argname, argval = argval)
        funccall += """
        _tpf = %(pfname)s
        _tobjclass = getattr(_tpf.h, '%(objtype)s')
        data_objects.append(_tobjclass(**_tobjargs))
        """ % dict(pfname = pfname, objtype = objtype)
        funccall = "\n".join((line.strip() for line in funccall.splitlines()))
        self.execute(funccall, hide = False)
        pf = self.locals['_tpf']

class ExtDirectParameterFileList(BottleDirectRouter):
    my_name = "ExtDirectParameterFileList"
    api_url = "pflist"

    def get_list_of_pfs(self):
        # Note that this instantiates the hierarchy.  This can be a costly
        # event.  However, we're going to assume that it's okay, if you have
        # decided to load up the parameter file.
        from yt.data_objects.static_output import _cached_pfs
        rv = []
        for fn, pf in sorted(_cached_pfs.items()):
            objs = []
            pf_varname = "_cached_pfs['%s']" % (fn)
            fields = []
            if pf._instantiated_hierarchy is not None: 
                fields = set(pf.h.field_list + pf.h.derived_field_list)
                fields = list(fields)
                fields.sort()
                for i,obj in enumerate(pf.h.objects):
                    try:
                        name = str(obj)
                    except ReferenceError:
                        continue
                    objs.append(dict(name=name, type=obj._type_name,
                                     filename = '', field_list = [],
                                     varname = "%s.h.objects[%s]" % (pf_varname, i)))
            rv.append( dict(name = str(pf), children = objs, filename=fn,
                            type = "parameter_file",
                            varname = pf_varname, field_list = fields) )
        return rv

def ext_load_script(filename):
    contents = open(filename).read()
    payload_handler = PayloadHandler()
    payload_handler.add_payload(
        {'type': 'cell_contents',
         'value': contents}
    )
    return

class PayloadLoggingHandler(logging.StreamHandler):
    def __init__(self, *args, **kwargs):
        logging.StreamHandler.__init__(self, *args, **kwargs)
        self.payload_handler = PayloadHandler()

    def emit(self, record):
        msg = self.format(record)
        self.payload_handler.add_payload(
            {'type':'logentry',
             'log_entry':msg})

if os.path.exists(os.path.expanduser("~/.yt/favicon.ico")):
    ico = os.path.expanduser("~/.yt/favicon.ico")
else:
    ico = os.path.join(local_dir, "html", "images", "favicon.ico")
@route("/favicon.ico", method="GET")
def _favicon_ico():
    response.headers['Content-Type'] = "image/x-icon"
    return open(ico).read()

class ExtProgressBar(object):
    def __init__(self, title, maxval):
        self.title = title
        self.maxval = maxval
        self.last = 0
        # Now we add a payload for the progress bar
        self.payload_handler = PayloadHandler()
        self.payload_handler.add_payload(
            {'type': 'widget',
             'widget_type': 'progressbar',
             'varname': None,
             'data': {'title':title}
            })

    def update(self, val):
        # An update is only meaningful if it's on the order of 1/100 or greater
        if ceil(100*self.last / self.maxval) + 1 == \
           floor(100*val / self.maxval) or val == self.maxval:
            self.last = val
            self.payload_handler.add_payload(
                {'type': 'widget_payload',
                 'widget_id': 'pbar_top',
                 'value': float(val) / self.maxval})

    def finish(self):
        self.payload_handler.add_payload(
            {'type': 'widget_payload',
             'widget_id': 'pbar_top',
             'value': -1})
