"""
Useful functions.  If non-original, see function for citation.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
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

import __builtin__
import time, types, signal, inspect, traceback, sys, pdb, os
import contextlib
import warnings, struct, subprocess
import numpy as np
from distutils import version
from math import floor, ceil

from yt.utilities.exceptions import *
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.definitions import inv_axis_names, axis_names, x_dict, y_dict
import yt.utilities.progressbar as pb
import yt.utilities.rpdb as rpdb
from collections import defaultdict
from functools import wraps

# Some functions for handling sequences and other types

def iterable(obj):
    """
    Grabbed from Python Cookbook / matploblib.cbook.  Returns true/false for
    *obj* iterable.
    """
    try: len(obj)
    except: return False
    return True

def ensure_list(obj):
    """
    This function ensures that *obj* is a list.  Typically used to convert a
    string to a list, for instance ensuring the *fields* as an argument is a
    list.
    """
    if obj is None:
        return [obj]
    if not isinstance(obj, types.ListType):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that *obj* is a numpy array. Typically used to
    convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        return obj
    elif isinstance(obj, (types.ListType, types.TupleType)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])

def ensure_tuple(obj):
    """
    This function ensures that *obj* is a tuple.  Typically used to convert
    scalar, list, or array arguments specified by a user in a context where
    we assume a tuple internally
    """
    if isintance(obj, types.TupleType):
        return obj
    elif isinstance(obj, (types.ListType, np.ndarray)):
        return tuple(obj)
    else:
        return (obj,)

def read_struct(f, fmt):
    """
    This reads a struct, and only that struct, from an open file.
    """
    s = f.read(struct.calcsize(fmt))
    return struct.unpack(fmt, s)

def just_one(obj):
    # If we have an iterable, sometimes we only want one item
    if hasattr(obj,'flat'):
        return obj.flat[0]
    elif iterable(obj):
        return obj[0]
    return obj

# Taken from
# http://www.goldb.org/goldblog/2008/02/06/PythonConvertSecsIntoHumanReadableTimeStringHHMMSS.aspx
def humanize_time(secs):
    """
    Takes *secs* and returns a nicely formatted string
    """
    mins, secs = divmod(secs, 60)
    hours, mins = divmod(mins, 60)
    return '%02d:%02d:%02d' % (hours, mins, secs)

#
# Some function wrappers that come in handy once in a while
#

# we use the resource module to get the memory page size

try:
    import resource
except ImportError:
    pass

def get_memory_usage():
    """
    Returning resident size in megabytes
    """
    pid = os.getpid()
    try:
        pagesize = resource.getpagesize()
    except NameError:
        return -1024
    status_file = "/proc/%s/statm" % (pid)
    if not os.path.isfile(status_file):
        return -1024
    line = open(status_file).read()
    size, resident, share, text, library, data, dt = [int(i) for i in line.split()]
    return resident * pagesize / (1024 * 1024) # return in megs

def time_execution(func):
    r"""
    Decorator for seeing how long a given function takes, depending on whether
    or not the global 'yt.timefunctions' config parameter is set.
    """
    @wraps(func)
    def wrapper(*arg, **kw):
        t1 = time.time()
        res = func(*arg, **kw)
        t2 = time.time()
        mylog.debug('%s took %0.3f s', func.func_name, (t2-t1))
        return res
    from yt.config import ytcfg
    if ytcfg.getboolean("yt","timefunctions") == True:
        return wrapper
    else:
        return func

def print_tb(func):
    """
    This function is used as a decorate on a function to have the calling stack
    printed whenever that function is entered.

    This can be used like so:

    .. code-block:: python

       @print_tb
       def some_deeply_nested_function(...):

    """
    @wraps(func)
    def run_func(*args, **kwargs):
        traceback.print_stack()
        return func(*args, **kwargs)
    return run_func

def rootonly(func):
    """
    This is a decorator that, when used, will only call the function on the
    root processor and then broadcast the results of the function to all other
    processors.

    This can be used like so:

    .. code-block:: python

       @rootonly
       def some_root_only_function(...):

    """
    from yt.config import ytcfg
    @wraps(func)
    def check_parallel_rank(*args, **kwargs):
        if ytcfg.getint("yt","__topcomm_parallel_rank") > 0:
            return
        return func(*args, **kwargs)
    return check_parallel_rank

def rootloginfo(*args):
    from yt.config import ytcfg
    if ytcfg.getint("yt", "__topcomm_parallel_rank") > 0: return
    mylog.info(*args)

def deprecate(func):
    """
    This decorator issues a deprecation warning.

    This can be used like so:

    .. code-block:: python

       @deprecate
       def some_really_old_function(...):

    """
    @wraps(func)
    def run_func(*args, **kwargs):
        warnings.warn("%s has been deprecated and may be removed without notice!" \
                % func.func_name, DeprecationWarning, stacklevel=2)
        func(*args, **kwargs)
    return run_func

def pdb_run(func):
    """
    This decorator inserts a pdb session on top of the call-stack into a
    function.

    This can be used like so:

    .. code-block:: python

       @pdb_run
       def some_function_to_debug(...):

    """
    @wraps(func)
    def wrapper(*args, **kw):
        pdb.runcall(func, *args, **kw)
    return wrapper

__header = """
== Welcome to the embedded IPython Shell ==

   You are currently inside the function:
     %(fname)s

   Defined in:
     %(filename)s:%(lineno)s
"""

def insert_ipython(num_up=1):
    """
    Placed inside a function, this will insert an IPython interpreter at that
    current location.  This will enabled detailed inspection of the current
    exeuction environment, as well as (optional) modification of that environment.
    *num_up* refers to how many frames of the stack get stripped off, and
    defaults to 1 so that this function itself is stripped off.
    """

    import IPython
    if version.LooseVersion(IPython.__version__) <= version.LooseVersion('0.10'):
        api_version = '0.10'
    else:
        api_version = '0.11'

    stack = inspect.stack()
    frame = inspect.stack()[num_up]
    loc = frame[0].f_locals.copy()
    glo = frame[0].f_globals
    dd = dict(fname = frame[3], filename = frame[1],
              lineno = frame[2])
    if api_version == '0.10':
        ipshell = IPython.Shell.IPShellEmbed()
        ipshell(header = __header % dd,
                local_ns = loc, global_ns = glo)
    else:
        from IPython.config.loader import Config
        cfg = Config()
        cfg.InteractiveShellEmbed.local_ns = loc
        cfg.InteractiveShellEmbed.global_ns = glo
        IPython.embed(config=cfg, banner2 = __header % dd)
        from IPython.frontend.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed(config=cfg)

    del ipshell


#
# Our progress bar types and how to get one
#

class DummyProgressBar(object):
    # This progressbar gets handed if we don't
    # want ANY output
    def __init__(self, *args, **kwargs):
        return
    def update(self, *args, **kwargs):
        return
    def finish(self, *args, **kwargs):
        return

class ParallelProgressBar(object):
    # This is just a simple progress bar
    # that prints on start/stop
    def __init__(self, title, maxval):
        self.title = title
        mylog.info("Starting '%s'", title)
    def update(self, *args, **kwargs):
        return
    def finish(self):
        mylog.info("Finishing '%s'", self.title)

class GUIProgressBar(object):
    def __init__(self, title, maxval):
        import wx
        self.maxval = maxval
        self.last = 0
        self._pbar = wx.ProgressDialog("Working...",
                    title, maximum=maxval,
                    style=wx.PD_REMAINING_TIME|wx.PD_ELAPSED_TIME|wx.PD_APP_MODAL)
    def update(self, val):
        # An update is only meaningful if it's on the order of 1/100 or greater
        if ceil(100*self.last / self.maxval) + 1 == \
           floor(100*val / self.maxval) or val == self.maxval:
            self._pbar.Update(val)
            self.last = val
    def finish(self):
        self._pbar.Destroy()

def get_pbar(title, maxval):
    """
    This returns a progressbar of the most appropriate type, given a *title*
    and a *maxval*.
    """
    maxval = max(maxval, 1)
    from yt.config import ytcfg
    if ytcfg.getboolean("yt", "suppressStreamLogging") or \
       "__IPYTHON__" in dir(__builtin__) or \
       ytcfg.getboolean("yt", "__withintesting"):
        return DummyProgressBar()
    elif ytcfg.getboolean("yt", "__withinreason"):
        from yt.gui.reason.extdirect_repl import ExtProgressBar
        return ExtProgressBar(title, maxval)
    elif ytcfg.getboolean("yt", "__parallel"):
        return ParallelProgressBar(title, maxval)
    widgets = [ title,
            pb.Percentage(), ' ',
            pb.Bar(marker=pb.RotatingMarker()),
            ' ', pb.ETA(), ' ']
    pbar = pb.ProgressBar(widgets=widgets,
                          maxval=maxval).start()
    return pbar

def only_on_root(func, *args, **kwargs):
    """
    This function accepts a *func*, a set of *args* and *kwargs* and then only
    on the root processor calls the function.  All other processors get "None"
    handed back.
    """
    from yt.config import ytcfg
    if kwargs.pop("global_rootonly", False):
        cfg_option = "__global_parallel_rank"
    else:
        cfg_option = "__topcomm_parallel_rank"
    if not ytcfg.getboolean("yt","__parallel"):
        return func(*args,**kwargs)
    if ytcfg.getint("yt", cfg_option) > 0: return
    return func(*args, **kwargs)

#
# Our signal and traceback handling functions
#

def signal_print_traceback(signo, frame):
    print traceback.print_stack(frame)

def signal_problem(signo, frame):
    raise RuntimeError()

def signal_ipython(signo, frame):
    insert_ipython(2)

def paste_traceback(exc_type, exc, tb):
    """
    This is a traceback handler that knows how to paste to the pastebin.
    Should only be used in sys.excepthook.
    """
    sys.__excepthook__(exc_type, exc, tb)
    import xmlrpclib, cStringIO
    p = xmlrpclib.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
    s = cStringIO.StringIO()
    traceback.print_exception(exc_type, exc, tb, file=s)
    s = s.getvalue()
    ret = p.pastes.newPaste('pytb', s, None, '', '', True)
    print
    print "Traceback pasted to http://paste.yt-project.org/show/%s" % (ret)
    print

def paste_traceback_detailed(exc_type, exc, tb):
    """
    This is a traceback handler that knows how to paste to the pastebin.
    Should only be used in sys.excepthook.
    """
    import xmlrpclib, cStringIO, cgitb
    s = cStringIO.StringIO()
    handler = cgitb.Hook(format="text", file = s)
    handler(exc_type, exc, tb)
    s = s.getvalue()
    print s
    p = xmlrpclib.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
    ret = p.pastes.newPaste('text', s, None, '', '', True)
    print
    print "Traceback pasted to http://paste.yt-project.org/show/%s" % (ret)
    print

def traceback_writer_hook(file_suffix = ""):
    def write_to_file(exc_type, exc, tb):
        sys.__excepthook__(exc_type, exc, tb)
        fn = "yt_traceback%s" % file_suffix
        f = open(fn, "w")
        traceback.print_exception(exc_type, exc, tb, file=f)
        print "Wrote traceback to %s" % fn
    return write_to_file

_ss = "fURbBUUBE0cLXgETJnZgJRMXVhVGUQpQAUBuehQMUhJWRFFRAV1ERAtBXw1dAxMLXT4zXBFfABNN\nC0ZEXw1YUURHCxMXVlFERwxWCQw=\n"
def _rdbeta(key):
    import itertools, base64
    enc_s = base64.decodestring(_ss)
    dec_s = ''.join([ chr(ord(a) ^ ord(b)) for a, b in zip(enc_s, itertools.cycle(key)) ])
    print dec_s

#
# Some exceptions
#

class NoCUDAException(Exception):
    pass

class YTEmptyClass(object):
    pass

def update_hg(path, skip_rebuild = False):
    from mercurial import hg, ui, commands
    f = open(os.path.join(path, "yt_updater.log"), "a")
    u = ui.ui()
    u.pushbuffer()
    config_fn = os.path.join(path, ".hg", "hgrc")
    print "Reading configuration from ", config_fn
    u.readconfig(config_fn)
    repo = hg.repository(u, path)
    commands.pull(u, repo)
    f.write(u.popbuffer())
    f.write("\n\n")
    u.pushbuffer()
    commands.identify(u, repo)
    if "+" in u.popbuffer():
        print "Can't rebuild modules by myself."
        print "You will have to do this yourself.  Here's a sample commands:"
        print
        print "    $ cd %s" % (path)
        print "    $ hg up"
        print "    $ %s setup.py develop" % (sys.executable)
        return 1
    print "Updating the repository"
    f.write("Updating the repository\n\n")
    commands.update(u, repo, check=True)
    if skip_rebuild: return
    f.write("Rebuilding modules\n\n")
    p = subprocess.Popen([sys.executable, "setup.py", "build_ext", "-i"], cwd=path,
                        stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    stdout, stderr = p.communicate()
    f.write(stdout)
    f.write("\n\n")
    if p.returncode:
        print "BROKEN: See %s" % (os.path.join(path, "yt_updater.log"))
        sys.exit(1)
    f.write("Successful!\n")
    print "Updated successfully."

def get_hg_version(path):
    from mercurial import hg, ui, commands
    u = ui.ui()
    u.pushbuffer()
    repo = hg.repository(u, path)
    commands.identify(u, repo)
    return u.popbuffer()

def get_yt_version():
    try:
        from yt.__hg_version__ import hg_version
        return hg_version
    except ImportError:
        pass
    import pkg_resources
    yt_provider = pkg_resources.get_provider("yt")
    path = os.path.dirname(yt_provider.module_path)
    version = get_hg_version(path)[:12]
    return version

def get_version_stack():
    import numpy, matplotlib, h5py
    version_info = {}
    version_info['yt'] = get_yt_version()
    version_info['numpy'] = numpy.version.version
    version_info['matplotlib'] = matplotlib.__version__
    version_info['h5py'] = h5py.version.version
    return version_info

def get_script_contents():
    stack = inspect.stack()
    top_frame = inspect.stack()[-1]
    finfo = inspect.getframeinfo(top_frame[0])
    if finfo[2] != "<module>": return None
    if not os.path.exists(finfo[0]): return None
    try:
        contents = open(finfo[0]).read()
    except:
        contents = None
    return contents

# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass = True):
    import urllib, urllib2
    uri = 'https://api.bitbucket.org/1.0/%s/' % endpoint
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.urlencode(data)
    req = urllib2.Request(uri, data)
    if use_pass:
        username = raw_input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = '%s:%s' % (username, password)
        req.add_header('Authorization', 'Basic %s' % base64.b64encode(upw).strip())
    return urllib2.urlopen(req).read()

def get_yt_supp():
    supp_path = os.path.join(os.environ["YT_DEST"], "src",
                             "yt-supplemental")
    # Now we check that the supplemental repository is checked out.
    if not os.path.isdir(supp_path):
        print
        print "*** The yt-supplemental repository is not checked ***"
        print "*** out.  I can do this for you, but because this ***"
        print "*** is a delicate act, I require you to respond   ***"
        print "*** to the prompt with the word 'yes'.            ***"
        print
        response = raw_input("Do you want me to try to check it out? ")
        if response != "yes":
            print
            print "Okay, I understand.  You can check it out yourself."
            print "This command will do it:"
            print
            print "$ hg clone http://hg.yt-project.org/yt-supplemental/ ",
            print "%s" % (supp_path)
            print
            sys.exit(1)
        rv = commands.clone(uu,
                "http://hg.yt-project.org/yt-supplemental/", supp_path)
        if rv:
            print "Something has gone wrong.  Quitting."
            sys.exit(1)
    # Now we think we have our supplemental repository.
    return supp_path

def fix_length(length, pf):
    if isinstance(length, (list, tuple)) and len(length) == 2 and \
       isinstance(length[1], types.StringTypes):
       length = length[0]/pf[length[1]]
    return length

@contextlib.contextmanager
def parallel_profile(prefix):
    import cProfile
    from yt.config import ytcfg
    fn = "%s_%04i_%04i.cprof" % (prefix,
                ytcfg.getint("yt", "__topcomm_parallel_size"),
                ytcfg.getint("yt", "__topcomm_parallel_rank"))
    p = cProfile.Profile()
    p.enable()
    yield fn
    p.disable()
    p.dump_stats(fn)

def get_num_threads():
    from .config import ytcfg
    nt = ytcfg.getint("yt","numthreads")
    if nt < 0:
        return os.environ.get("OMP_NUM_THREADS", 0)
    return nt

def fix_axis(axis):
    return inv_axis_names.get(axis, axis)

def get_image_suffix(name):
    suffix = os.path.splitext(name)[1]
    return suffix if suffix in ['.png', '.eps', '.ps', '.pdf'] else ''
