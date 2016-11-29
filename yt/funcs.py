"""
Useful functions.  If non-original, see function for citation.



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import errno
from yt.extern.six import string_types
from yt.extern.six.moves import input
import time
import inspect
import traceback
import sys
import pdb
import os
import re
import contextlib
import warnings
import struct
import subprocess
import numpy as np
import itertools
import base64
import numpy
import matplotlib
import getpass
from distutils.version import LooseVersion
from math import floor, ceil
from numbers import Number as numeric_type

from yt.extern.six.moves import urllib
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.exceptions import YTInvalidWidthError
from yt.extern.tqdm import tqdm
from yt.units.yt_array import YTArray, YTQuantity
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
    if not isinstance(obj, list):
        return [obj]
    return obj

def ensure_numpy_array(obj):
    """
    This function ensures that *obj* is a numpy array. Typically used to
    convert scalar, list or tuple argument passed to functions using Cython.
    """
    if isinstance(obj, np.ndarray):
        if obj.shape == ():
            return np.array([obj])
        # We cast to ndarray to catch ndarray subclasses
        return np.array(obj)
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    else:
        return np.asarray([obj])

def ensure_tuple(obj):
    """
    This function ensures that *obj* is a tuple.  Typically used to convert
    scalar, list, or array arguments specified by a user in a context where
    we assume a tuple internally
    """
    if isinstance(obj, tuple):
        return obj
    elif isinstance(obj, (list, np.ndarray)):
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
        if isinstance(obj, YTArray):
            return YTQuantity(obj.flat[0], obj.units, registry=obj.units.registry)
        return obj.flat[0]
    elif iterable(obj):
        return obj[0]
    return obj


def compare_dicts(dict1, dict2):
    if not set(dict1) <= set(dict2):
        return False
    for key in dict1.keys():
        if dict1[key] is not None and dict2[key] is not None:
            if isinstance(dict1[key], dict):
                if compare_dicts(dict1[key], dict2[key]):
                    continue
                else:
                    return False
            try:
                comparison = (dict1[key] == dict2[key]).all()
            except AttributeError:
                comparison = (dict1[key] == dict2[key])
            if not comparison:
                return False
    return True

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

def get_memory_usage(subtract_share = False):
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
    if subtract_share: resident -= share
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
        mylog.debug('%s took %0.3f s', func.__name__, (t2-t1))
        return res
    from yt.config import ytcfg
    if ytcfg.getboolean("yt","timefunctions") is True:
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

def deprecate(replacement):
    def real_deprecate(func):
        """
        This decorator issues a deprecation warning.

        This can be used like so:

        .. code-block:: python

        @deprecate("new_function")
        def some_really_old_function(...):

        """
        @wraps(func)
        def run_func(*args, **kwargs):
            message = "%s has been deprecated and may be removed without notice!"
            if replacement is not None:
                message += " Use %s instead." % replacement
            warnings.warn(message % func.__name__, DeprecationWarning,
                          stacklevel=2)
            func(*args, **kwargs)
        return run_func
    return real_deprecate

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

def get_ipython_api_version():
    import IPython
    if LooseVersion(IPython.__version__) <= LooseVersion('0.10'):
        api_version = '0.10'
    elif LooseVersion(IPython.__version__) <= LooseVersion('1.0'):
        api_version = '0.11'
    else:
        api_version = '1.0'

    return api_version

def insert_ipython(num_up=1):
    """
    Placed inside a function, this will insert an IPython interpreter at that
    current location.  This will enabled detailed inspection of the current
    exeuction environment, as well as (optional) modification of that environment.
    *num_up* refers to how many frames of the stack get stripped off, and
    defaults to 1 so that this function itself is stripped off.
    """
    import IPython
    api_version = get_ipython_api_version()

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
        try:
            from traitlets.config.loader import Config
        except ImportError:
            from IPython.config.loader import Config
        cfg = Config()
        cfg.InteractiveShellEmbed.local_ns = loc
        cfg.InteractiveShellEmbed.global_ns = glo
        IPython.embed(config=cfg, banner2 = __header % dd)
        if api_version == '0.11':
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
        else:
            from IPython.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed(config=cfg)

    del ipshell


#
# Our progress bar types and how to get one
#

class TqdmProgressBar(object):
    # This is a drop in replacement for pbar
    # called tqdm
    def __init__(self,title, maxval):
        self._pbar = tqdm(leave=True, total=maxval, desc=title)
        self.i = 0
    def update(self, i=None):
        if i is None:
            i = self.i + 1
        n = i - self.i
        self.i = i
        self._pbar.update(n)
    def finish(self):
        self._pbar.close()

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

def get_pbar(title, maxval, parallel=False):
    """
    This returns a progressbar of the most appropriate type, given a *title*
    and a *maxval*.
    """
    maxval = max(maxval, 1)
    from yt.config import ytcfg
    if ytcfg.getboolean("yt", "suppressStreamLogging") or \
       ytcfg.getboolean("yt", "__withintesting"):
        return DummyProgressBar()
    elif ytcfg.getboolean("yt", "__parallel"):
        # If parallel is True, update progress on root only.
        if parallel:
            if is_root():
                return TqdmProgressBar(title, maxval)
            else:
                return DummyProgressBar()
        else:
            return ParallelProgressBar(title, maxval)
    pbar = TqdmProgressBar(title,maxval)
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

def is_root():
    """
    This function returns True if it is on the root processor of the
    topcomm and False otherwise.
    """
    from yt.config import ytcfg
    cfg_option = "__topcomm_parallel_rank"
    if not ytcfg.getboolean("yt","__parallel"):
        return True
    if ytcfg.getint("yt", cfg_option) > 0:
        return False
    return True


#
# Our signal and traceback handling functions
#

def signal_print_traceback(signo, frame):
    print(traceback.print_stack(frame))

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
    from yt.extern.six.moves import StringIO, xmlrpc_client
    p = xmlrpc_client.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
    s = StringIO()
    traceback.print_exception(exc_type, exc, tb, file=s)
    s = s.getvalue()
    ret = p.pastes.newPaste('pytb', s, None, '', '', True)
    print()
    print("Traceback pasted to http://paste.yt-project.org/show/%s" % (ret))
    print()

def paste_traceback_detailed(exc_type, exc, tb):
    """
    This is a traceback handler that knows how to paste to the pastebin.
    Should only be used in sys.excepthook.
    """
    import cgitb
    from yt.extern.six.moves import StringIO, xmlrpc_client
    s = StringIO()
    handler = cgitb.Hook(format="text", file = s)
    handler(exc_type, exc, tb)
    s = s.getvalue()
    print(s)
    p = xmlrpc_client.ServerProxy(
            "http://paste.yt-project.org/xmlrpc/",
            allow_none=True)
    ret = p.pastes.newPaste('text', s, None, '', '', True)
    print()
    print("Traceback pasted to http://paste.yt-project.org/show/%s" % (ret))
    print()

_ss = "fURbBUUBE0cLXgETJnZgJRMXVhVGUQpQAUBuehQMUhJWRFFRAV1ERAtBXw1dAxMLXT4zXBFfABNN\nC0ZEXw1YUURHCxMXVlFERwxWCQw=\n"
def _rdbeta(key):
    enc_s = base64.decodestring(_ss)
    dec_s = ''.join([ chr(ord(a) ^ ord(b)) for a, b in zip(enc_s, itertools.cycle(key)) ])
    print(dec_s)

#
# Some exceptions
#

class NoCUDAException(Exception):
    pass

class YTEmptyClass(object):
    pass

def update_hg(path, skip_rebuild = False):
    try:
        import hglib
    except ImportError:
        print("Updating requires python-hglib to be installed.")
        print("Try: pip install python-hglib")
        return -1
    f = open(os.path.join(path, "yt_updater.log"), "a")
    with hglib.open(path) as repo:
        repo.pull()
        ident = repo.identify().decode("utf-8")
        if "+" in ident:
            print("Can't rebuild modules by myself.")
            print("You will have to do this yourself.  Here's a sample commands:")
            print("")
            print("    $ cd %s" % (path))
            print("    $ hg up")
            print("    $ %s setup.py develop" % (sys.executable))
            return 1
        print("Updating the repository")
        f.write("Updating the repository\n\n")
        repo.update(check=True)
        f.write("Updated from %s to %s\n\n" % (ident, repo.identify()))
        if skip_rebuild: return
        f.write("Rebuilding modules\n\n")
        p = subprocess.Popen([sys.executable, "setup.py", "build_ext", "-i"],
                             cwd=path, stdout = subprocess.PIPE,
                             stderr = subprocess.STDOUT)
        stdout, stderr = p.communicate()
        f.write(stdout.decode('utf-8'))
        f.write("\n\n")
        if p.returncode:
            print("BROKEN: See %s" % (os.path.join(path, "yt_updater.log")))
            sys.exit(1)
        f.write("Successful!\n")
        print("Updated successfully.")

def get_hg_version(path):
    try:
        import hglib
    except ImportError:
        print("Updating and precise version information requires ")
        print("python-hglib to be installed.")
        print("Try: pip install python-hglib")
        return None
    try:
        with hglib.open(path) as repo:
            return repo.identify()
    except hglib.error.ServerError:
        # path is not an hg repository
        return None

def get_yt_version():
    try:
        from yt.__hg_version__ import hg_version
        return hg_version
    except ImportError:
        pass
    import pkg_resources
    yt_provider = pkg_resources.get_provider("yt")
    path = os.path.dirname(yt_provider.module_path)
    version = get_hg_version(path)
    if version is None:
        return version
    else:
        return version[:12].strip().decode('utf-8')

def get_version_stack():
    version_info = {}
    version_info['yt'] = get_yt_version()
    version_info['numpy'] = numpy.version.version
    version_info['matplotlib'] = matplotlib.__version__
    return version_info

def get_script_contents():
    top_frame = inspect.stack()[-1]
    finfo = inspect.getframeinfo(top_frame[0])
    if finfo[2] != "<module>": return None
    if not os.path.exists(finfo[0]): return None
    try:
        contents = open(finfo[0]).read()
    except:
        contents = None
    return contents

def download_file(url, filename):
    class MyURLopener(urllib.request.FancyURLopener):
        def http_error_default(self, url, fp, errcode, errmsg, headers):
            raise RuntimeError("Attempt to download file from %s failed with error %s: %s." % \
              (url, errcode, errmsg))
    fn, h = MyURLopener().retrieve(url, filename)
    return fn

# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass = True):
    uri = 'https://api.bitbucket.org/1.0/%s/' % endpoint
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.parse.urlencode(data)
    req = urllib.request.Request(uri, data)
    if use_pass:
        username = input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = '%s:%s' % (username, password)
        req.add_header('Authorization', 'Basic %s' % base64.b64encode(upw).strip())
    return urllib.request.urlopen(req).read()

def get_yt_supp():
    import hglib
    supp_path = os.path.join(os.environ["YT_DEST"], "src",
                             "yt-supplemental")
    # Now we check that the supplemental repository is checked out.
    if not os.path.isdir(supp_path):
        print()
        print("*** The yt-supplemental repository is not checked ***")
        print("*** out.  I can do this for you, but because this ***")
        print("*** is a delicate act, I require you to respond   ***")
        print("*** to the prompt with the word 'yes'.            ***")
        print()
        response = input("Do you want me to try to check it out? ")
        if response != "yes":
            print()
            print("Okay, I understand.  You can check it out yourself.")
            print("This command will do it:")
            print()
            print("$ hg clone http://bitbucket.org/yt_analysis/yt-supplemental/ ", end=' ')
            print("%s" % (supp_path))
            print()
            sys.exit(1)
        rv = hglib.clone("http://bitbucket.org/yt_analysis/yt-supplemental/", 
                         supp_path)
        if rv:
            print("Something has gone wrong.  Quitting.")
            sys.exit(1)
    # Now we think we have our supplemental repository.
    return supp_path

def fix_length(length, ds=None):
    assert ds is not None
    if ds is not None:
        registry = ds.unit_registry
    else:
        registry = None
    if isinstance(length, YTArray):
        if registry is not None:
            length.units.registry = registry
        return length.in_units("code_length")
    if isinstance(length, numeric_type):
        return YTArray(length, 'code_length', registry=registry)
    length_valid_tuple = isinstance(length, (list, tuple)) and len(length) == 2
    unit_is_string = isinstance(length[1], string_types)
    if length_valid_tuple and unit_is_string:
        return YTArray(*length, registry=registry)
    else:
        raise RuntimeError("Length %s is invalid" % str(length))

@contextlib.contextmanager
def parallel_profile(prefix):
    r"""A context manager for profiling parallel code execution using cProfile

    This is a simple context manager that automatically profiles the execution
    of a snippet of code.

    Parameters
    ----------
    prefix : string
        A string name to prefix outputs with.

    Examples
    --------

    >>> with parallel_profile('my_profile'):
    ...     yt.PhasePlot(ds.all_data(), 'density', 'temperature', 'cell_mass')
    """
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

def fix_axis(axis, ds):
    return ds.coordinates.axis_id.get(axis, axis)

def get_image_suffix(name):
    suffix = os.path.splitext(name)[1]
    return suffix if suffix in ['.png', '.eps', '.ps', '.pdf'] else ''

def get_output_filename(name, keyword, suffix):
    r"""Return an appropriate filename for output.

    With a name provided by the user, this will decide how to 
    appropriately name the output file by the following rules:
    1. if name is None, the filename will be the keyword plus 
       the suffix.
    2. if name ends with "/", assume name is a directory and 
       the file will be named name/(keyword+suffix).  If the
       directory does not exist, first try to create it and
       raise an exception if an error occurs.
    3. if name does not end in the suffix, add the suffix.
    
    Parameters
    ----------
    name : str
        A filename given by the user.
    keyword : str
        A default filename prefix if name is None.
    suffix : str
        Suffix that must appear at end of the filename.
        This will be added if not present.

    Examples
    --------

    >>> print get_output_filename(None, "Projection_x", ".png")
    Projection_x.png
    >>> print get_output_filename("my_file", "Projection_x", ".png")
    my_file.png
    >>> print get_output_filename("my_file/", "Projection_x", ".png")
    my_file/Projection_x.png
    
    """
    if name is None:
        name = keyword
    name = os.path.expanduser(name)
    if name[-1] == os.sep and not os.path.isdir(name):
        try:
            os.mkdir(name)
        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise
    if os.path.isdir(name):
        name = os.path.join(name, keyword)
    if not name.endswith(suffix):
        name += suffix
    return name

def ensure_dir_exists(path):
    r"""Create all directories in path recursively in a parallel safe manner"""
    my_dir = os.path.dirname(path)
    if not my_dir:
        return
    if not os.path.exists(my_dir):
        only_on_root(os.makedirs, my_dir)

def ensure_dir(path):
    r"""Parallel safe directory maker."""
    if not os.path.exists(path):
        only_on_root(os.makedirs, path)
    return path

def validate_width_tuple(width):
    if not iterable(width) or len(width) != 2:
        raise YTInvalidWidthError("width (%s) is not a two element tuple" % width)
    if not isinstance(width[0], numeric_type) and isinstance(width[1], string_types):
        msg = "width (%s) is invalid. " % str(width)
        msg += "Valid widths look like this: (12, 'au')"
        raise YTInvalidWidthError(msg)

def camelcase_to_underscore(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

def set_intersection(some_list):
    if len(some_list) == 0: return set([])
    # This accepts a list of iterables, which we get the intersection of.
    s = set(some_list[0])
    for l in some_list[1:]:
        s.intersection_update(l)
    return s

@contextlib.contextmanager
def memory_checker(interval = 15, dest = None):
    r"""This is a context manager that monitors memory usage.

    Parameters
    ----------
    interval : int
        The number of seconds between printing the current memory usage in
        gigabytes of the current Python interpreter.

    Examples
    --------

    >>> with memory_checker(10):
    ...     arr = np.zeros(1024*1024*1024, dtype="float64")
    ...     time.sleep(15)
    ...     del arr
    """
    import threading
    if dest is None:
        dest = sys.stdout
    class MemoryChecker(threading.Thread):
        def __init__(self, event, interval):
            self.event = event
            self.interval = interval
            threading.Thread.__init__(self)

        def run(self):
            while not self.event.wait(self.interval):
                print("MEMORY: %0.3e gb" % (get_memory_usage()/1024.), file=dest)

    e = threading.Event()
    mem_check = MemoryChecker(e, interval)
    mem_check.start()
    try:
        yield
    finally:
        e.set()


def deprecated_class(cls):
    @wraps(cls)
    def _func(*args, **kwargs):
        # Note we use SyntaxWarning because by default, DeprecationWarning is
        # not shown.
        warnings.warn(
            "This usage is deprecated.  Please use %s instead." % cls.__name__,
            SyntaxWarning, stacklevel=2)
        return cls(*args, **kwargs)
    return _func

def enable_plugins():
    """Forces the plugins file to be parsed.

    This plugin file is a means of creating custom fields, quantities,
    data objects, colormaps, and other code classes and objects to be used
    in yt scripts without modifying the yt source directly.

    The file must be located at ``$HOME/.yt/my_plugins.py``.

    Warning: when you use this function, your script will only be reproducible
    if you also provide the ``my_plugins.py`` file.
    """
    import yt
    from yt.fields.my_plugin_fields import my_plugins_fields
    from yt.config import ytcfg
    my_plugin_name = ytcfg.get("yt","pluginfilename")
    # We assume that it is with respect to the $HOME/.yt directory
    if os.path.isfile(my_plugin_name):
        _fn = my_plugin_name
    else:
        _fn = os.path.expanduser("~/.yt/%s" % my_plugin_name)
    if os.path.isfile(_fn):
        mylog.info("Loading plugins from %s", _fn)
        execdict = yt.__dict__.copy()
        execdict['add_field'] = my_plugins_fields.add_field
        with open(_fn) as f:
            code = compile(f.read(), _fn, 'exec')
            exec(code, execdict)

def fix_unitary(u):
    if u == '1':
        return 'unitary'
    else:
        return u

def get_hash(infile, algorithm='md5', BLOCKSIZE=65536):
    """Generate file hash without reading in the entire file at once.

    Original code licensed under MIT.  Source:
    http://pythoncentral.io/hashing-files-with-python/

    Parameters
    ----------
    infile : str
        File of interest (including the path).
    algorithm : str (optional)
        Hash algorithm of choice. Defaults to 'md5'.
    BLOCKSIZE : int (optional)
        How much data in bytes to read in at once.

    Returns
    -------
    hash : str
        The hash of the file.

    Examples
    --------
    >>> import yt.funcs as funcs
    >>> funcs.get_hash('/path/to/test.png')
    'd38da04859093d430fa4084fd605de60'

    """
    import hashlib

    try:
        hasher = getattr(hashlib, algorithm)()
    except:
        raise NotImplementedError("'%s' not available!  Available algorithms: %s" %
                                  (algorithm, hashlib.algorithms))

    filesize   = os.path.getsize(infile)
    iterations = int(float(filesize)/float(BLOCKSIZE))

    pbar = get_pbar('Generating %s hash' % algorithm, iterations)

    iter = 0
    with open(infile,'rb') as f:
        buf = f.read(BLOCKSIZE)
        while len(buf) > 0:
            hasher.update(buf)
            buf = f.read(BLOCKSIZE)
            iter += 1
            pbar.update(iter)
        pbar.finish()

    return hasher.hexdigest()

def get_brewer_cmap(cmap):
    """Returns a colorbrewer colormap from palettable"""
    try:
        import brewer2mpl
    except ImportError:
        brewer2mpl = None
    try:
        import palettable
    except ImportError:
        palettable = None
    if palettable is not None:
        bmap = palettable.colorbrewer.get_map(*cmap)
    elif brewer2mpl is not None:
        warnings.warn("Using brewer2mpl colormaps is deprecated. "
                      "Please install the successor to brewer2mpl, "
                      "palettable, with `pip install palettable`. "
                      "Colormap tuple names remain unchanged.")
        bmap = brewer2mpl.get_map(*cmap)
    else:
        raise RuntimeError(
            "Please install palettable to use colorbrewer colormaps")
    return bmap.get_mpl_colormap(N=cmap[2])

def get_requests():
    try:
        import requests
    except ImportError:
        requests = None
    return requests

@contextlib.contextmanager
def dummy_context_manager(*args, **kwargs):
    yield

def matplotlib_style_context(style_name=None, after_reset=True):
    """Returns a context manager for controlling matplotlib style.

    Arguments are passed to matplotlib.style.context() if specified. Defaults
    to setting "classic" style, after resetting to the default config parameters.

    On older matplotlib versions (<=1.5.0) where matplotlib.style isn't
    available, returns a dummy context manager.
    """
    if style_name is None:
        style_name = 'classic'
    try:
        import matplotlib.style
        if style_name in matplotlib.style.available:
            return matplotlib.style.context(style_name, after_reset=after_reset)
    except ImportError:
        pass
    return dummy_context_manager()

def setdefaultattr(obj, name, value):
    """Set attribute with *name* on *obj* with *value* if it doesn't exist yet

    Analogous to dict.setdefault
    """
    if not hasattr(obj, name):
        setattr(obj, name, value)
    return getattr(obj, name)
