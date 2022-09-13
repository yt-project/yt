import base64
import builtins
import contextlib
import copy
import errno
import getpass
import glob
import inspect
import itertools
import os
import pdb
import re
import struct
import subprocess
import sys
import time
import traceback
import urllib.parse
import urllib.request
from collections import UserDict
from functools import lru_cache, wraps
from numbers import Number as numeric_type
from typing import Any, Callable, Type

import numpy as np
from more_itertools import always_iterable, collapse, first
from packaging.version import Version
from tqdm import tqdm

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.units import YTArray, YTQuantity
from yt.utilities.exceptions import YTInvalidWidthError
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _requests as requests

# Some functions for handling sequences and other types


def is_sequence(obj):
    """
    Grabbed from Python Cookbook / matplotlib.cbook.  Returns true/false for

    Parameters
    ----------
    obj : iterable
    """
    try:
        len(obj)
        return True
    except TypeError:
        return False


def iter_fields(field_or_fields):
    """
    Create an iterator for field names, specified as single strings or tuples(fname,
    ftype) alike.
    This can safely be used in places where we accept a single field or a list as input.

    Parameters
    ----------
    field_or_fields: str, tuple(str, str), or any iterable of the previous types.

    Examples
    --------

    >>> fields = ("gas", "density")
    >>> for field in iter_fields(fields):
    ...     print(field)
    density

    >>> fields = ("gas", "density")
    >>> for field in iter_fields(fields):
    ...     print(field)
    ('gas', 'density')

    >>> fields = [("gas", "density"), ("gas", "temperature"), ("index", "dx")]
    >>> for field in iter_fields(fields):
    ...     print(field)
    density
    temperature
    ('index', 'dx')
    """
    return always_iterable(field_or_fields, base_type=(tuple, str, bytes))


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


def read_struct(f, fmt):
    """
    This reads a struct, and only that struct, from an open file.
    """
    s = f.read(struct.calcsize(fmt))
    return struct.unpack(fmt, s)


def just_one(obj):
    # If we have an iterable, sometimes we only want one item
    return first(collapse(obj))


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
                comparison = np.array_equal(dict1[key], dict2[key])
            except TypeError:
                comparison = dict1[key] == dict2[key]
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
    return "%02d:%02d:%02d" % (hours, mins, secs)


#
# Some function wrappers that come in handy once in a while
#

# we use the resource module to get the memory page size

try:
    import resource
except ImportError:
    pass


def get_memory_usage(subtract_share=False):
    """
    Returning resident size in megabytes
    """
    pid = os.getpid()
    try:
        pagesize = resource.getpagesize()
    except NameError:
        return -1024
    status_file = f"/proc/{pid}/statm"
    if not os.path.isfile(status_file):
        return -1024
    with open(status_file) as fh:
        line = fh.read()
    size, resident, share, text, library, data, dt = (int(i) for i in line.split())
    if subtract_share:
        resident -= share
    return resident * pagesize / (1024 * 1024)  # return in megs


def time_execution(func):
    r"""
    Decorator for seeing how long a given function takes, depending on whether
    or not the global 'yt.time_functions' config parameter is set.
    """

    @wraps(func)
    def wrapper(*arg, **kw):
        t1 = time.time()
        res = func(*arg, **kw)
        t2 = time.time()
        mylog.debug("%s took %0.3f s", func.__name__, (t2 - t1))
        return res

    if ytcfg.get("yt", "time_functions"):
        return wrapper
    else:
        return func


def print_tb(func):
    """
    This function is used as a decorate on a function to have the calling stack
    printed whenever that function is entered.

    This can be used like so:

    >>> @print_tb
    ... def some_deeply_nested_function(*args, **kwargs):
    ...     ...

    """

    @wraps(func)
    def run_func(*args, **kwargs):
        traceback.print_stack()
        return func(*args, **kwargs)

    return run_func


def rootonly(func):
    """
    This is a decorator that, when used, will only call the function on the
    root processor.

    This can be used like so:

    .. code-block:: python

       @rootonly
       def some_root_only_function(*args, **kwargs):
           ...
    """

    @wraps(func)
    def check_parallel_rank(*args, **kwargs):
        if ytcfg.get("yt", "internals", "topcomm_parallel_rank") > 0:
            return
        return func(*args, **kwargs)

    return check_parallel_rank


def pdb_run(func):
    """
    This decorator inserts a pdb session on top of the call-stack into a
    function.

    This can be used like so:

    >>> @pdb_run
    ... def some_function_to_debug(*args, **kwargs):
    ...     ...

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
    execution environment, as well as (optional) modification of that environment.
    *num_up* refers to how many frames of the stack get stripped off, and
    defaults to 1 so that this function itself is stripped off.
    """
    import IPython
    from IPython.terminal.embed import InteractiveShellEmbed

    try:
        from traitlets.config.loader import Config
    except ImportError:
        from IPython.config.loader import Config

    frame = inspect.stack()[num_up]
    loc = frame[0].f_locals.copy()
    glo = frame[0].f_globals
    dd = dict(fname=frame[3], filename=frame[1], lineno=frame[2])
    cfg = Config()
    cfg.InteractiveShellEmbed.local_ns = loc
    cfg.InteractiveShellEmbed.global_ns = glo
    IPython.embed(config=cfg, banner2=__header % dd)
    ipshell = InteractiveShellEmbed(config=cfg)

    del ipshell


#
# Our progress bar types and how to get one
#


class TqdmProgressBar:
    # This is a drop in replacement for pbar
    # called tqdm
    def __init__(self, title, maxval):
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


class DummyProgressBar:
    # This progressbar gets handed if we don't
    # want ANY output
    def __init__(self, *args, **kwargs):
        return

    def update(self, *args, **kwargs):
        return

    def finish(self, *args, **kwargs):
        return


def get_pbar(title, maxval):
    """
    This returns a progressbar of the most appropriate type, given a *title*
    and a *maxval*.
    """
    maxval = max(maxval, 1)

    if (
        ytcfg.get("yt", "suppress_stream_logging")
        or ytcfg.get("yt", "internals", "within_testing")
        or maxval == 1
        or not is_root()
    ):
        return DummyProgressBar()
    return TqdmProgressBar(title, maxval)


def only_on_root(func, *args, **kwargs):
    """
    This function accepts a *func*, a set of *args* and *kwargs* and then only
    on the root processor calls the function.  All other processors get "None"
    handed back.
    """
    if kwargs.pop("global_rootonly", False):
        cfg_option = "global_parallel_rank"
    else:
        cfg_option = "topcomm_parallel_rank"
    if not ytcfg.get("yt", "internals", "parallel"):
        return func(*args, **kwargs)
    if ytcfg.get("yt", "internals", cfg_option) > 0:
        return
    return func(*args, **kwargs)


def is_root():
    """
    This function returns True if it is on the root processor of the
    topcomm and False otherwise.
    """
    if not ytcfg.get("yt", "internals", "parallel"):
        return True
    return ytcfg.get("yt", "internals", "topcomm_parallel_rank") == 0


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
    import xmlrpc.client
    from io import StringIO

    p = xmlrpc.client.ServerProxy(
        "http://paste.yt-project.org/xmlrpc/", allow_none=True
    )
    s = StringIO()
    traceback.print_exception(exc_type, exc, tb, file=s)
    s = s.getvalue()
    ret = p.pastes.newPaste("pytb", s, None, "", "", True)
    print()
    print(f"Traceback pasted to http://paste.yt-project.org/show/{ret}")
    print()


def paste_traceback_detailed(exc_type, exc, tb):
    """
    This is a traceback handler that knows how to paste to the pastebin.
    Should only be used in sys.excepthook.
    """
    import cgitb
    import xmlrpc.client
    from io import StringIO

    s = StringIO()
    handler = cgitb.Hook(format="text", file=s)
    handler(exc_type, exc, tb)
    s = s.getvalue()
    print(s)
    p = xmlrpc.client.ServerProxy(
        "http://paste.yt-project.org/xmlrpc/", allow_none=True
    )
    ret = p.pastes.newPaste("text", s, None, "", "", True)
    print()
    print(f"Traceback pasted to http://paste.yt-project.org/show/{ret}")
    print()


_ss = "fURbBUUBE0cLXgETJnZgJRMXVhVGUQpQAUBuehQMUhJWRFFRAV1ERAtBXw1dAxMLXT4zXBFfABNN\nC0ZEXw1YUURHCxMXVlFERwxWCQw=\n"


def _rdbeta(key):
    enc_s = base64.decodestring(_ss)
    dec_s = "".join(chr(ord(a) ^ ord(b)) for a, b in zip(enc_s, itertools.cycle(key)))
    print(dec_s)


#
# Some exceptions
#


class NoCUDAException(Exception):
    pass


class YTEmptyClass:
    pass


def update_git(path):
    try:
        import git
    except ImportError:
        print("Updating and precise version information requires ")
        print("gitpython to be installed.")
        print("Try: python -m pip install gitpython")
        return -1
    with open(os.path.join(path, "yt_updater.log"), "a") as f:
        repo = git.Repo(path)
        if repo.is_dirty(untracked_files=True):
            print("Changes have been made to the yt source code so I won't ")
            print("update the code. You will have to do this yourself.")
            print("Here's a set of sample commands:")
            print("")
            print(f"    $ cd {path}")
            print("    $ git stash")
            print("    $ git checkout main")
            print("    $ git pull")
            print("    $ git stash pop")
            print(f"    $ {sys.executable} setup.py develop")
            print("")
            return 1
        if repo.active_branch.name != "main":
            print("yt repository is not tracking the main branch so I won't ")
            print("update the code. You will have to do this yourself.")
            print("Here's a set of sample commands:")
            print("")
            print(f"    $ cd {path}")
            print("    $ git checkout main")
            print("    $ git pull")
            print(f"    $ {sys.executable} setup.py develop")
            print("")
            return 1
        print("Updating the repository")
        f.write("Updating the repository\n\n")
        old_version = repo.git.rev_parse("HEAD", short=12)
        try:
            remote = repo.remotes.yt_upstream
        except AttributeError:
            remote = repo.create_remote(
                "yt_upstream", url="https://github.com/yt-project/yt"
            )
            remote.fetch()
        main = repo.heads.main
        main.set_tracking_branch(remote.refs.main)
        main.checkout()
        remote.pull()
        new_version = repo.git.rev_parse("HEAD", short=12)
        f.write(f"Updated from {old_version} to {new_version}\n\n")
        rebuild_modules(path, f)
    print("Updated successfully")


def rebuild_modules(path, f):
    f.write("Rebuilding modules\n\n")
    p = subprocess.Popen(
        [sys.executable, "setup.py", "build_ext", "-i"],
        cwd=path,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    stdout, stderr = p.communicate()
    f.write(stdout.decode("utf-8"))
    f.write("\n\n")
    if p.returncode:
        print(f"BROKEN: See {os.path.join(path, 'yt_updater.log')}")
        sys.exit(1)
    f.write("Successful!\n")


def get_git_version(path):
    try:
        import git
    except ImportError:
        print("Updating and precise version information requires ")
        print("gitpython to be installed.")
        print("Try: python -m pip install gitpython")
        return None
    try:
        repo = git.Repo(path)
        return repo.git.rev_parse("HEAD", short=12)
    except git.InvalidGitRepositoryError:
        # path is not a git repository
        return None


def get_yt_version():
    import pkg_resources

    yt_provider = pkg_resources.get_provider("yt")
    path = os.path.dirname(yt_provider.module_path)
    version = get_git_version(path)
    if version is None:
        return version
    else:
        v_str = version[:12].strip()
        if hasattr(v_str, "decode"):
            v_str = v_str.decode("utf-8")
        return v_str


def get_version_stack():
    import matplotlib

    version_info = {}
    version_info["yt"] = get_yt_version()
    version_info["numpy"] = np.version.version
    version_info["matplotlib"] = matplotlib.__version__
    return version_info


def get_script_contents():
    top_frame = inspect.stack()[-1]
    finfo = inspect.getframeinfo(top_frame[0])
    if finfo[2] != "<module>":
        return None
    if not os.path.exists(finfo[0]):
        return None
    try:
        contents = open(finfo[0]).read()
    except Exception:
        contents = None
    return contents


def download_file(url, filename):
    try:
        return fancy_download_file(url, filename, requests)
    except ImportError:
        # fancy_download_file requires requests
        return simple_download_file(url, filename)


def fancy_download_file(url, filename, requests=None):
    response = requests.get(url, stream=True)
    total_length = response.headers.get("content-length")

    with open(filename, "wb") as fh:
        if total_length is None:
            fh.write(response.content)
        else:
            blocksize = 4 * 1024**2
            iterations = int(float(total_length) / float(blocksize))

            pbar = get_pbar(
                "Downloading %s to %s " % os.path.split(filename)[::-1], iterations
            )
            iteration = 0
            for chunk in response.iter_content(chunk_size=blocksize):
                fh.write(chunk)
                iteration += 1
                pbar.update(iteration)
            pbar.finish()
    return filename


def simple_download_file(url, filename):
    class MyURLopener(urllib.request.FancyURLopener):
        def http_error_default(self, url, fp, errcode, errmsg, headers):
            raise RuntimeError(
                "Attempt to download file from %s failed with error %s: %s."
                % (url, errcode, errmsg)
            )

    fn, h = MyURLopener().retrieve(url, filename)
    return fn


# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass=True):
    uri = f"https://api.bitbucket.org/1.0/{endpoint}/"
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.parse.urlencode(data)
    req = urllib.request.Request(uri, data)
    if use_pass:
        username = input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = f"{username}:{password}"
        req.add_header("Authorization", f"Basic {base64.b64encode(upw).strip()}")
    return urllib.request.urlopen(req).read()


def fix_length(length, ds):
    registry = ds.unit_registry
    if isinstance(length, YTArray):
        if registry is not None:
            length.units.registry = registry
        return length.in_units("code_length")
    if isinstance(length, numeric_type):
        return YTArray(length, "code_length", registry=registry)
    length_valid_tuple = isinstance(length, (list, tuple)) and len(length) == 2
    unit_is_string = isinstance(length[1], str)
    length_is_number = isinstance(length[0], numeric_type) and not isinstance(
        length[0], YTArray
    )
    if length_valid_tuple and unit_is_string and length_is_number:
        return YTArray(*length, registry=registry)
    else:
        raise RuntimeError(f"Length {str(length)} is invalid")


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

    >>> from yt import PhasePlot
    >>> from yt.testing import fake_random_ds
    >>> fields = ("density", "temperature", "cell_mass")
    >>> units = ("g/cm**3", "K", "g")
    >>> ds = fake_random_ds(16, fields=fields, units=units)
    >>> with parallel_profile("my_profile"):
    ...     plot = PhasePlot(ds.all_data(), *fields)
    """
    import cProfile

    fn = "%s_%04i_%04i.cprof" % (
        prefix,
        ytcfg.get("yt", "internals", "topcomm_parallel_size"),
        ytcfg.get("yt", "internals", "topcomm_parallel_rank"),
    )
    p = cProfile.Profile()
    p.enable()
    yield fn
    p.disable()
    p.dump_stats(fn)


def get_num_threads():
    from .config import ytcfg

    nt = ytcfg.get("yt", "num_threads")
    if nt < 0:
        return os.environ.get("OMP_NUM_THREADS", 0)
    return nt


def fix_axis(axis, ds):
    return ds.coordinates.axis_id.get(axis, axis)


def get_output_filename(name, keyword, suffix):
    r"""Return an appropriate filename for output.

    With a name provided by the user, this will decide how to appropriately name the
    output file by the following rules:

    1. if name is None, the filename will be the keyword plus the suffix.
    2. if name ends with "/" (resp "\" on Windows), assume name is a directory and the
       file will be named name/(keyword+suffix).  If the directory does not exist, first
       try to create it and raise an exception if an error occurs.
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

    >>> get_output_filename(None, "Projection_x", ".png")
    'Projection_x.png'
    >>> get_output_filename("my_file", "Projection_x", ".png")
    'my_file.png'
    >>> get_output_filename("my_dir/", "Projection_x", ".png")
    'my_dir/Projection_x.png'

    """
    if name is None:
        name = keyword
    name = os.path.expanduser(name)
    if name.endswith(os.sep) and not os.path.isdir(name):
        ensure_dir(name)
    if os.path.isdir(name):
        name = os.path.join(name, keyword)
    if not name.endswith(suffix):
        name += suffix
    return name


def ensure_dir_exists(path):
    r"""Create all directories in path recursively in a parallel safe manner"""
    my_dir = os.path.dirname(path)
    # If path is a file in the current directory, like "test.txt", then my_dir
    # would be an empty string, resulting in FileNotFoundError when passed to
    # ensure_dir. Let's avoid that.
    if my_dir:
        ensure_dir(my_dir)


def ensure_dir(path):
    r"""Parallel safe directory maker."""
    if os.path.exists(path):
        return path

    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    return path


def validate_width_tuple(width):
    if not is_sequence(width) or len(width) != 2:
        raise YTInvalidWidthError(f"width ({width}) is not a two element tuple")
    is_numeric = isinstance(width[0], numeric_type)
    length_has_units = isinstance(width[0], YTArray)
    unit_is_string = isinstance(width[1], str)
    if not is_numeric or length_has_units and unit_is_string:
        msg = f"width ({str(width)}) is invalid. "
        msg += "Valid widths look like this: (12, 'au')"
        raise YTInvalidWidthError(msg)


_first_cap_re = re.compile("(.)([A-Z][a-z]+)")
_all_cap_re = re.compile("([a-z0-9])([A-Z])")


@lru_cache(maxsize=128, typed=False)
def camelcase_to_underscore(name):
    s1 = _first_cap_re.sub(r"\1_\2", name)
    return _all_cap_re.sub(r"\1_\2", s1).lower()


def set_intersection(some_list):
    if len(some_list) == 0:
        return set()
    # This accepts a list of iterables, which we get the intersection of.
    s = set(some_list[0])
    for l in some_list[1:]:
        s.intersection_update(l)
    return s


@contextlib.contextmanager
def memory_checker(interval=15, dest=None):
    r"""This is a context manager that monitors memory usage.

    Parameters
    ----------
    interval : int
        The number of seconds between printing the current memory usage in
        gigabytes of the current Python interpreter.

    Examples
    --------

    >>> with memory_checker(10):
    ...     arr = np.zeros(1024 * 1024 * 1024, dtype="float64")
    ...     time.sleep(15)
    ...     del arr
    MEMORY: -1.000e+00 gb
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
                print(f"MEMORY: {get_memory_usage() / 1024.0:0.3e} gb", file=dest)

    e = threading.Event()
    mem_check = MemoryChecker(e, interval)
    mem_check.start()
    try:
        yield
    finally:
        e.set()


def enable_plugins(plugin_filename=None):
    """Forces a plugin file to be parsed.

    A plugin file is a means of creating custom fields, quantities,
    data objects, colormaps, and other code classes and objects to be used
    in yt scripts without modifying the yt source directly.

    If ``plugin_filename`` is omitted, this function will look for a plugin file at
    ``$HOME/.config/yt/my_plugins.py``, which is the preferred behaviour for a
    system-level configuration.

    Warning: a script using this function will only be reproducible if your plugin
    file is shared with it.
    """
    import yt
    from yt.config import config_dir, ytcfg
    from yt.fields.my_plugin_fields import my_plugins_fields

    if plugin_filename is not None:
        _fn = plugin_filename
        if not os.path.isfile(_fn):
            raise FileNotFoundError(_fn)
    else:
        # Determine global plugin location. By decreasing priority order:
        # - absolute path
        # - CONFIG_DIR
        # - obsolete config dir.
        my_plugin_name = ytcfg.get("yt", "plugin_filename")
        for base_prefix in ("", config_dir()):
            if os.path.isfile(os.path.join(base_prefix, my_plugin_name)):
                _fn = os.path.join(base_prefix, my_plugin_name)
                break
        else:
            raise FileNotFoundError("Could not find a global system plugin file.")

    mylog.info("Loading plugins from %s", _fn)
    ytdict = yt.__dict__
    execdict = ytdict.copy()
    execdict["add_field"] = my_plugins_fields.add_field
    with open(_fn) as f:
        code = compile(f.read(), _fn, "exec")
        exec(code, execdict, execdict)
    ytnamespace = list(ytdict.keys())
    for k in execdict.keys():
        if k not in ytnamespace:
            if callable(execdict[k]):
                setattr(yt, k, execdict[k])


def subchunk_count(n_total, chunk_size):
    handled = 0
    while handled < n_total:
        tr = min(n_total - handled, chunk_size)
        yield tr
        handled += tr


def fix_unitary(u):
    if u == "1":
        return "unitary"
    else:
        return u


def get_hash(infile, algorithm="md5", BLOCKSIZE=65536):
    """Generate file hash without reading in the entire file at once.

    Original code licensed under MIT.  Source:
    https://www.pythoncentral.io/hashing-files-with-python/

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
    >>> from tempfile import NamedTemporaryFile
    >>> with NamedTemporaryFile() as file:
    ...     get_hash(file.name)
    'd41d8cd98f00b204e9800998ecf8427e'
    """
    import hashlib

    try:
        hasher = getattr(hashlib, algorithm)()
    except AttributeError as e:
        raise NotImplementedError(
            f"'{algorithm}' not available!  Available algorithms: {hashlib.algorithms}"
        ) from e

    filesize = os.path.getsize(infile)
    iterations = int(float(filesize) / float(BLOCKSIZE))

    pbar = get_pbar(f"Generating {algorithm} hash", iterations)

    iter = 0
    with open(infile, "rb") as f:
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
        issue_deprecation_warning(
            "Using brewer2mpl colormaps is deprecated. "
            "Please install the successor to brewer2mpl, "
            "palettable, with `pip install palettable`. "
            "Colormap tuple names remain unchanged.",
            since="3.3",
            removal="4.2",
        )
        bmap = brewer2mpl.get_map(*cmap)
    else:
        raise RuntimeError("Please install palettable to use colorbrewer colormaps")
    return bmap.get_mpl_colormap(N=cmap[2])


@contextlib.contextmanager
def dummy_context_manager(*args, **kwargs):
    yield


def matplotlib_style_context(style_name=None, after_reset=False):
    """Returns a context manager for controlling matplotlib style.

    Arguments are passed to matplotlib.style.context() if specified. Defaults
    to setting "classic" style, after resetting to the default config parameters.

    On older matplotlib versions (<=1.5.0) where matplotlib.style isn't
    available, returns a dummy context manager.
    """
    if style_name is None:
        import matplotlib

        style_name = {"mathtext.fontset": "cm"}
        if Version(matplotlib.__version__) >= Version("3.3.0"):
            style_name["mathtext.fallback"] = "cm"
        else:
            style_name["mathtext.fallback_to_cm"] = True
    try:
        import matplotlib.style

        return matplotlib.style.context(style_name, after_reset=after_reset)
    except ImportError:
        pass
    return dummy_context_manager()


interactivity = False

"""Sets the condition that interactive backends can be used."""


def toggle_interactivity():
    global interactivity
    interactivity = not interactivity
    if interactivity:
        if "__IPYTHON__" in dir(builtins):
            import IPython

            shell = IPython.get_ipython()
            shell.magic("matplotlib")
        else:
            import matplotlib

            matplotlib.interactive(True)


def get_interactivity():
    return interactivity


def setdefaultattr(obj, name, value):
    """Set attribute with *name* on *obj* with *value* if it doesn't exist yet

    Analogous to dict.setdefault
    """
    if not hasattr(obj, name):
        setattr(obj, name, value)
    return getattr(obj, name)


def parse_h5_attr(f, attr):
    """A Python3-safe function for getting hdf5 attributes.

    If an attribute is supposed to be a string, this will return it as such.
    """
    val = f.attrs.get(attr, None)
    if isinstance(val, bytes):
        return val.decode("utf8")
    else:
        return val


def obj_length(v):
    if is_sequence(v):
        return len(v)
    else:
        # If something isn't iterable, we return 0
        # to signify zero length (aka a scalar).
        return 0


def array_like_field(data, x, field):
    field = data._determine_fields(field)[0]
    if isinstance(field, tuple):
        finfo = data.ds._get_field_info(field[0], field[1])
    else:
        finfo = data.ds._get_field_info(field)
    if finfo.sampling_type == "particle":
        units = finfo.output_units
    else:
        units = finfo.units
    if isinstance(x, YTArray):
        arr = copy.deepcopy(x)
        arr.convert_to_units(units)
        return arr
    if isinstance(x, np.ndarray):
        return data.ds.arr(x, units)
    else:
        return data.ds.quan(x, units)


def validate_3d_array(obj):
    if not is_sequence(obj) or len(obj) != 3:
        raise TypeError(
            "Expected an array of size (3,), received '%s' of "
            "length %s" % (str(type(obj)).split("'")[1], len(obj))
        )


def validate_float(obj):
    """Validates if the passed argument is a float value.

    Raises an exception if `obj` is a single float value
    or a YTQuantity of size 1.

    Parameters
    ----------
    obj : Any
        Any argument which needs to be checked for a single float value.

    Raises
    ------
    TypeError
        Raised if `obj` is not a single float value or YTQunatity

    Examples
    --------
    >>> validate_float(1)
    >>> validate_float(1.50)
    >>> validate_float(YTQuantity(1, "cm"))
    >>> validate_float((1, "cm"))
    >>> validate_float([1, 1, 1])
    Traceback (most recent call last):
    ...
    TypeError: Expected a numeric value (or size-1 array), received 'list' of length 3

    >>> validate_float([YTQuantity(1, "cm"), YTQuantity(2, "cm")])
    Traceback (most recent call last):
    ...
    TypeError: Expected a numeric value (or size-1 array), received 'list' of length 2
    """
    if isinstance(obj, tuple):
        if (
            len(obj) != 2
            or not isinstance(obj[0], numeric_type)
            or not isinstance(obj[1], str)
        ):
            raise TypeError(
                "Expected a numeric value (or tuple of format "
                "(float, String)), received an inconsistent tuple "
                "'%s'." % str(obj)
            )
        else:
            return
    if is_sequence(obj) and (len(obj) != 1 or not isinstance(obj[0], numeric_type)):
        raise TypeError(
            "Expected a numeric value (or size-1 array), "
            "received '%s' of length %s" % (str(type(obj)).split("'")[1], len(obj))
        )


def validate_sequence(obj):
    if obj is not None and not is_sequence(obj):
        raise TypeError(
            "Expected an iterable object,"
            " received '%s'" % str(type(obj)).split("'")[1]
        )


def validate_field_key(key):
    if (
        isinstance(key, tuple)
        and len(key) == 2
        and all(isinstance(_, str) for _ in key)
    ):
        return
    raise TypeError(
        "Expected a 2-tuple of strings formatted as\n"
        "(field or particle type, field name)\n"
        f"Received invalid field key: {key}, with type {type(key)}"
    )


def validate_object(obj, data_type):
    if obj is not None and not isinstance(obj, data_type):
        raise TypeError(
            "Expected an object of '%s' type, received '%s'"
            % (str(data_type).split("'")[1], str(type(obj)).split("'")[1])
        )


def validate_axis(ds, axis):
    if ds is not None:
        valid_axis = ds.coordinates.axis_name.keys()
    else:
        valid_axis = [0, 1, 2, "x", "y", "z", "X", "Y", "Z"]
    if axis not in valid_axis:
        raise TypeError(
            "Expected axis of int or char type (can be %s), "
            "received '%s'." % (list(valid_axis), axis)
        )


def validate_center(center):
    if isinstance(center, str):
        c = center.lower()
        if (
            c not in ["c", "center", "m", "max", "min"]
            and not c.startswith("max_")
            and not c.startswith("min_")
        ):
            raise TypeError(
                "Expected 'center' to be in ['c', 'center', "
                "'m', 'max', 'min'] or the prefix to be "
                "'max_'/'min_', received '%s'." % center
            )
    elif not isinstance(center, (numeric_type, YTQuantity)) and not is_sequence(center):
        raise TypeError(
            "Expected 'center' to be a numeric object of type "
            "list/tuple/np.ndarray/YTArray/YTQuantity, "
            "received '%s'." % str(type(center)).split("'")[1]
        )


def sglob(pattern):
    """
    Return the results of a glob through the sorted() function.
    """
    return sorted(glob.glob(pattern))


def dictWithFactory(factory: Callable[[Any], Any]) -> Type:
    """
    Create a dictionary class with a default factory function.
    Contrary to `collections.defaultdict`, the factory takes
    the missing key as input parameter.

    Parameters
    ----------
    factory : callable(key) -> value
        The factory to call when hitting a missing key

    Returns
    -------
    DictWithFactory class
        A class to create new dictionaries handling missing keys.
    """

    issue_deprecation_warning(
        "yt.funcs.dictWithFactory will be removed in a future version of yt, please do not rely on it. "
        "If you need it, copy paste this function from yt's source code",
        since="4.1",
    )

    class DictWithFactory(UserDict):
        def __init__(self, *args, **kwargs):
            self.factory = factory
            super().__init__(*args, **kwargs)

        def __missing__(self, key):
            val = self.factory(key)
            self[key] = val
            return val

    return DictWithFactory


def levenshtein_distance(seq1, seq2, max_dist=None):
    """
    Compute the levenshtein distance between seq1 and seq2.
    From https://stackabuse.com/levenshtein-distance-and-text-similarity-in-python/

    Parameters
    ----------
    seq1 : str
    seq2 : str
        The strings to compute the distance between
    max_dist : integer
        If not None, maximum distance returned (see notes).

    Returns
    -------
    The Levenshtein distance as an integer.

    Notes
    -----
    This computes the Levenshtein distance, i.e. the number of edits to change
    seq1 into seq2. If a maximum distance is passed, the algorithm will stop as soon
    as the number of edits goes above the value. This allows for an earlier break
    and speeds calculations up.
    """
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    if max_dist is None:
        max_dist = max(size_x, size_y)

    if abs(size_x - size_y) > max_dist:
        return max_dist + 1
    matrix = np.zeros((size_x, size_y), dtype=int)
    for x in range(size_x):
        matrix[x, 0] = x
    for y in range(size_y):
        matrix[0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x - 1] == seq2[y - 1]:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1, matrix[x - 1, y - 1], matrix[x, y - 1] + 1
                )
            else:
                matrix[x, y] = min(
                    matrix[x - 1, y] + 1, matrix[x - 1, y - 1] + 1, matrix[x, y - 1] + 1
                )

        # Early break: the minimum distance is already larger than
        # maximum allow value, can return safely.
        if matrix[x].min() > max_dist:
            return max_dist + 1
    return matrix[size_x - 1, size_y - 1]


def validate_moment(moment, weight_field):
    if moment == 2 and weight_field is None:
        raise ValueError(
            "Cannot compute the second moment of a projection if weight_field=None!"
        )
    if moment not in [1, 2]:
        raise ValueError(
            "Weighted projections can only be made of averages "
            "(moment = 1) or standard deviations (moment = 2)!"
        )
