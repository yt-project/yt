"""
Useful functions.  If non-original, see function for citation.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

import time, types, signal, inspect, traceback, sys, pdb, rpdb
import warnings, os
import progressbar as pb
from math import floor, ceil
from yt.logger import ytLogger as mylog

def signal_print_traceback(signo, frame):
    print traceback.print_stack(frame)

def signal_problem(signo, frame):
    raise RuntimeError()

try:
    signal.signal(signal.SIGUSR1, signal_print_traceback)
    mylog.debug("SIGUSR1 registered for traceback printing")
    signal.signal(signal.SIGUSR2, signal_problem)
    mylog.debug("SIGUSR2 registered for RuntimeError")
except ValueError:  # Not in main thread
    pass

def paste_traceback(exc_type, exc, tb):
    sys.__excepthook__(exc_type, exc, tb)
    import xmlrpclib, cStringIO
    p = xmlrpclib.ServerProxy(
            "http://paste.enzotools.org/xmlrpc/",
            allow_none=True)
    s = cStringIO.StringIO()
    traceback.print_exception(exc_type, exc, tb, file=s)
    s = s.getvalue()
    ret = p.pastes.newPaste('pytb', s, None, '', '', True)
    print
    print "Traceback pasted to http://paste.enzotools.org/show/%s" % (ret)
    print

if "--paste" in sys.argv:
    sys.excepthook = paste_traceback
if "--rpdb" in sys.argv:
    sys.excepthook = rpdb.rpdb_excepthook
    del sys.argv[sys.argv.index("--rpdb")]

def blank_wrapper(f):
    return lambda a: a

try:
    from functools import wraps
except ImportError:
    wraps = blank_wrapper

def iterable(obj):
    """
    Grabbed from Python Cookbook / matploblib.cbook
    """
    try: len(obj)
    except: return False
    return True

def ensure_list(obj):
    if obj == None:
        return [obj]
    if not isinstance(obj, types.ListType):
        return [obj]
    return obj

def time_execution(func):
    """
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

def pdb_run(func):
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
    from IPython.Shell import IPShellEmbed
    stack = inspect.stack()
    frame = inspect.stack()[num_up]
    loc = frame[0].f_locals.copy()
    glo = frame[0].f_globals
    dd = dict(fname = frame[3], filename = frame[1],
              lineno = frame[2])
    ipshell = IPShellEmbed()
    ipshell(header = __header % dd,
            local_ns = loc, global_ns = glo)
    del ipshell

class DummyProgressBar:
    def __init__(self, *args, **kwargs):
        return
    def update(self, *args, **kwargs):
        return
    def finish(sefl, *args, **kwargs):
        return

class GUIProgressBar:
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

def just_one(obj):
    if hasattr(obj,'flat'):
        return obj.flat[0]
    elif iterable(obj):
        return obj[0]
    return obj

def get_pbar(title, maxval):
    from yt.config import ytcfg
    if ytcfg.getboolean("yt","inGui"):
        if maxval > ytcfg.getint("reason","minpbar"): # Arbitrary number
            return GUIProgressBar(title, maxval)
        else:
            return DummyProgressBar()
    elif ytcfg.getboolean("yt","suppressStreamLogging"):
        return DummyProgressBar()
    elif "SAGE_ROOT" in os.environ:
        try:
            from sage.server.support import EMBEDDED_MODE
            if EMBEDDED_MODE: return DummyProgressBar()
        except:
            pass
    widgets = [ title,
            pb.Percentage(), ' ',
            pb.Bar(marker=pb.RotatingMarker()),
            ' ', pb.ETA(), ' ']
    pbar = pb.ProgressBar(widgets=widgets,
                          maxval=maxval).start()
    return pbar

# Taken from
# http://www.goldb.org/goldblog/2008/02/06/PythonConvertSecsIntoHumanReadableTimeStringHHMMSS.aspx
def humanize_time(secs):
    mins, secs = divmod(secs, 60)
    hours, mins = divmod(mins, 60)
    return '%02d:%02d:%02d' % (hours, mins, secs)

class __defaultdict(dict):
    def __init__(self, func):
        self.__func = func
        dict.__init__(self)
    def __getitem__(self, key):
        if not self.has_key(key):
            self.__setitem__(key, self.__func())
        return dict.__getitem__(self, key)

import traceback
def print_tb(func):
    @wraps(func)
    def run_func(*args, **kwargs):
        traceback.print_stack()
        return func(*args, **kwargs)
    return run_func

try:
    from collections import defaultdict
except ImportError:
    defaultdict = __defaultdict

def rootonly(func):
    @wraps(func)
    def donothing(*args, **kwargs):
        return
    from yt.config import ytcfg
    if ytcfg.getint("yt","__parallel_rank") > 0: return donothing
    return func

def only_on_root(func, *args, **kwargs):
    from yt.config import ytcfg
    if not ytcfg.getboolean("yt","__parallel"):
        return func(*args,**kwargs)
    if ytcfg.getint("yt","__parallel_rank") > 0: return
    return func(*args, **kwargs)

def deprecate(func):
    @wraps(func)
    def run_func(*args, **kwargs):
        warnings.warn("%s has been deprecated and may be removed without notice!" \
                % func.func_name, DeprecationWarning, stacklevel=2)
        func(*args, **kwargs)
    return run_func
