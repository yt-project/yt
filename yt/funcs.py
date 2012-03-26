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

import time, types, signal, inspect, traceback, sys, pdb, os
import warnings, struct
from math import floor, ceil

from yt.utilities.exceptions import *
from yt.utilities.logger import ytLogger as mylog
import yt.utilities.progressbar as pb
import yt.utilities.rpdb as rpdb

# Some compatibility functions.  In the long run, these *should* disappear as
# we move toward newer python versions.  Most were implemented to get things
# running on DataStar.

# If we're running on python2.4, we need a 'wraps' function
def blank_wrapper(f):
    return lambda a: a

try:
    from functools import wraps
except ImportError:
    wraps = blank_wrapper

# We need to ensure that we have a defaultdict implementation

class __defaultdict(dict):
    def __init__(self, func):
        self.__func = func
        dict.__init__(self)
    def __getitem__(self, key):
        if not self.has_key(key):
            self.__setitem__(key, self.__func())
        return dict.__getitem__(self, key)

try:
    from collections import defaultdict
except ImportError:
    defaultdict = __defaultdict

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
    if obj == None:
        return [obj]
    if not isinstance(obj, types.ListType):
        return [obj]
    return obj

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
        return float(os.popen('ps -o rss= -p %d' % pid).read()) / 1024
    status_file = "/proc/%s/statm" % (pid)
    if not os.path.isfile(status_file):
        return float(os.popen('ps -o rss= -p %d' % pid).read()) / 1024
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

       @rootonly
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

       @rootonly
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
    if IPython.__version__.startswith("0.10"):
       api_version = '0.10'
    elif IPython.__version__.startswith("0.11") or \
         IPython.__version__.startswith("0.12"):
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
    if ytcfg.getboolean("yt","suppressStreamLogging"):
        return DummyProgressBar()
    elif ytcfg.getboolean("yt", "__parallel"):
        return ParallelProgressBar(title, maxval)
    elif "SAGE_ROOT" in os.environ:
        try:
            from sage.server.support import EMBEDDED_MODE
            if EMBEDDED_MODE: return DummyProgressBar()
        except:
            pass
    elif "CODENODE" in os.environ:
        return DummyProgressBar()
    elif ytcfg.getboolean("yt", "__withinreason"):
        from yt.gui.reason.extdirect_repl import ExtProgressBar
        return ExtProgressBar(title, maxval)
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

# We use two signals, SIGUSR1 and SIGUSR2.  In a non-threaded environment,
# we set up handlers to process these by printing the current stack and to
# raise a RuntimeError.  The latter can be used, inside pdb, to catch an error
# and then examine the current stack.
try:
    signal.signal(signal.SIGUSR1, signal_print_traceback)
    mylog.debug("SIGUSR1 registered for traceback printing")
    signal.signal(signal.SIGUSR2, signal_ipython)
    mylog.debug("SIGUSR2 registered for IPython Insertion")
except ValueError:  # Not in main thread
    pass

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

# If we recognize one of the arguments on the command line as indicating a
# different mechanism for handling tracebacks, we attach one of those handlers
# and remove the argument from sys.argv.
#
# This fallback is for Paraview:
if not hasattr(sys, 'argv') or sys.argv is None: sys.argv = []
# Now, we check.
if "--paste" in sys.argv:
    sys.excepthook = paste_traceback
    del sys.argv[sys.argv.index("--paste")]
elif "--paste-detailed" in sys.argv:
    sys.excepthook = paste_traceback_detailed
    del sys.argv[sys.argv.index("--paste-detailed")]
elif "--detailed" in sys.argv:
    import cgitb; cgitb.enable(format="text")
    del sys.argv[sys.argv.index("--detailed")]
elif "--rpdb" in sys.argv:
    sys.excepthook = rpdb.rpdb_excepthook
    del sys.argv[sys.argv.index("--rpdb")]
elif "--detailed" in sys.argv:
    import cgitb; cgitb.enable(format="text")
    del sys.argv[sys.argv.index("--detailed")]

#
# Some exceptions
#

class NoCUDAException(Exception):
    pass

class YTEmptyClass(object):
    pass
