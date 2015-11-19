#####
#
# This code is included in yt, but was not originally written by yt.  For the
# original source, please see this discussion:
#
#   https://groups.google.com/forum/#!topic/mpi4py/DYNXufzfKPA
#
#   and the repository:
#
#   https://github.com/langton/MPI_Import
#
#
# Inclusion with yt does not imply any change in license of the original code.
# Our local modifications, to fix recursive imports and to use mpi4py, are
# under released under the original code license.
#
#####


# This code is derived from knee.py, which was included in the Python
# 2.6 distribution.
#
# The modifications to this code are copyright (c) 2011, Lawrence
# Livermore National Security, LLC. Produced at the Lawrence Livermore
# National Laboratory. Written by Tim Kadich and Asher Langton
# <langton2@llnl.gov>. Released as LLNL-CODE-522751 under the name
# SmartImport.py, version 1.0. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# - Redistributions of source code must retain the above copyright
# notice, this list of conditions and the disclaimer below.
#
# - Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the disclaimer (as noted below)
#   in the documentation and/or other materials provided with the
#   distribution.
#
# - Neither the name of the LLNS/LLNL nor the names of its contributors
#   may be used to endorse or promote products derived from this
#   software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE
# LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Additional BSD Notice
#
# 1. This notice is required to be provided under our contract with the
# U.S. Department of Energy (DOE). This work was produced at Lawrence
# Livermore National Laboratory under Contract No. DE-AC52-07NA27344
# with the DOE.
#
# 2. Neither the United States Government nor Lawrence Livermore
# National Security, LLC nor any of their employees, makes any warranty,
# express or implied, or assumes any liability or responsibility for the
# accuracy, completeness, or usefulness of any information, apparatus,
# product, or process disclosed, or represents that its use would not
# infringe privately-owned rights.
#
# 3. Also, reference herein to any specific commercial products,
# process, or services by trade name, trademark, manufacturer or
# otherwise does not necessarily constitute or imply its endorsement,
# recommendation, or favoring by the United States Government or
# Lawrence Livermore National Security, LLC. The views and opinions of
# authors expressed herein do not necessarily state or reflect those of
# the United States Government or Lawrence Livermore National Security,
# LLC, and shall not be used for advertising or product endorsement
# purposes.

"""MPI_Import defines an mpi-aware import hook. The standard use of
this module is as follows:

   from MPI_Import import mpi_import
   with mpi_import():
      import foo
      import bar

Within the with block, the standard import statement is replaced by an
MPI-aware import statement. The rank 0 process finds the location of
each module to import, broadcasts the location, then all of the
processes load that module.

One CRITICAL detail: any code inside the mpi_import block must be
executed exactly the same on all of the MPI ranks. For example,
consider this:

def foo():
   import mpi
   if mpi.rank == 0:
      bar = someFunction()
   bar = mpi.bcast(bar,root=0)

def someFunction():
   import os
   return os.name

If foo() is called during the import process, then things may go very
wrong. If the os module hasn't been loaded, then the rank 0 process
will find os and broadcast its location. Since there's no
corresponding bcast for rank > 0, the other processes will receive
that broadcast instead of the broadcast for bar, resulting in
undefined behavior. Similarly, if rank >0 process encounters an import
that rank 0 does not encounter, that process will either hang waiting
for the bcast, or it will receive an out-of-order bcast.

The import hook provides a way to test whether we're using this
importer, which can be used to disable rank-asymmetric behavior in a
module import:

from yt.extern.six.moves import builtins
hasattr(builtins.__import__,"mpi_import")

This evaluates to True only when we're in an mpi_import() context
manager.

There are some situations where rank-dependent code may be necessary.
One such example is pyMPI's synchronizeQueuedOutput function, which
tends to cause deadlocks when it is executed inside an mpi_imported
module. In that case, we provide a hook to execute a function after
the mpi_import hook has been replaced by the standard import hook.
Here is an example showing the use of this feature:

# encapsulate the rank-asymmetric code in a function
def f():
    if mpi.rank == 0:
        doOneThing()
    else:
        doSomethingElse()

# Either importer is None (standard import) or it's a reference to
# the mpi_import object that owns the current importer.
from yt.extern.six.moves import builtins
importer = getattr(builtins.__import__,"mpi_import",None)
if importer:
    importer.callAfterImport(f)
else:
    # If we're using the standard import, then we'll execute the
    # code in f immediately
    f()

WARNING: the callAfterImport feature is not intended for casual use.
Usually it will be sufficient (and preferable) to either remove the
rank-asymmetric code or explicitly move it outside of the 'with
mpi_import' block. callAfterImport is provided for the (hopefully
rare!) cases where this does not suffice.


Some implementation details:

-This code is based on knee.py, which is an example of a pure Python
 hierarchical import that was included with Python 2.6 distributions.

-Python PEP 302 defines another way to override import by using finder
 and loader objects, which behave similarly to the imp.find_module and
 imp.load_module functions in __import_module__ below. Unfortunately,
 the implementation of PEP 302 is such that the path for the module
 has already been found by the time that the "finder" object is
 constructed, so it's not suitable for our purposes.

-This module uses pyMPI. It was originally designed with mpi4py, and
 switching back to mpi4py requires only minor modifications. To
 quickly substitute mpi4py for pyMPI, the 'import mpi' line below can
 be replaced with the following wrapper:

from mpi4py import MPI
class mpi(object):
    rank = MPI.COMM_WORLD.Get_rank()
    @staticmethod
    def bcast(obj=None,root=0):
        return MPI.COMM_WORLD.bcast(obj,root)

-An alternate version of this module had rank 0 perform all of the
 lookups, and then broadcast the locations all-at-once when that
 process reached the end of the context manager. This was somewhat
 faster than the current implementation, but was prone to deadlock
 when loading modules containing MPI synchronization points.

-The 'level' parameter to the import hook is not handled correctly; we
 treat it as if it were -1 (try relative and absolute imports). For
 more information about the level parameter, run 'help(__import__)'.
"""
from __future__ import print_function

import sys
import imp
import types
from yt.extern.six.moves import builtins
from mpi4py import MPI


class mpi(object):
    rank = MPI.COMM_WORLD.Get_rank()
    @staticmethod
    def bcast(obj=None,root=0):
        return MPI.COMM_WORLD.bcast(obj,root)

class mpi_import(object):
    def __enter__(self):
        imp.acquire_lock()
        __import_hook__.mpi_import = self
        self.__funcs = []
        self.original_import = builtins.__import__
        builtins.__import__ = __import_hook__

    def __exit__(self,type,value,traceback):
        builtins.__import__ = self.original_import
        __import_hook__.mpi_import = None
        imp.release_lock()
        for f in self.__funcs:
            f()

    def callAfterImport(self,f):
        "Add f to the list of functions to call on exit"
        if not isinstance(f, types.FunctionType):
            raise TypeError("Argument must be a function!")
        self.__funcs.append(f)


# The remaining code is for internal use only. Do not explicitly call
# call any of the following functions.

# Replacement for __import__(). Taken from knee.py; unmodified except for the
# (unused) level parameter.
def __import_hook__(name, globals=None, locals=None, fromlist=None, level=-1):
    # TODO: handle level parameter correctly. For now, we'll ignore
    # it and try both absolute and relative imports.
    parent = __determine_parent__(globals, level)
    q, tail = __find_head_package__(parent, name)
    m = __load_tail__(q, tail)
    if not fromlist:
        return q
    if hasattr(m, "__path__"):
        __ensure_fromlist__(m, fromlist)
    return m

# __import_module__ is the only part of knee.py with non-trivial changes.
# The MPI rank 0 process handles the lookup and broadcasts the location to
# the others. This must be called synchronously, at least in the case that
# 'fqname' is not already in sys.modules.
def __import_module__(partname, fqname, parent):
    fqname = fqname.rstrip(".")
    try:
        return sys.modules[fqname]
    except KeyError:
        pass
    fp = None         # module's file
    pathname = None   # module's location
    stuff = None      # tuple of (suffix,mode,type) for the module
    ierror = False    # are we propagating an import error from rank 0?

    # Start with the lookup on rank 0. The other processes will be waiting
    # on a broadcast, so we need to send one even if we're bailing out due
    # to an import error.
    if mpi.rank == 0:
        try:
            fp, pathname, stuff = imp.find_module(partname,
                                                  parent and parent.__path__)
        except ImportError:
            ierror = True
            return None
        finally:
            pathname,stuff,ierror = mpi.bcast((pathname,stuff,ierror))
    else:
        pathname,stuff,ierror = mpi.bcast((pathname,stuff,ierror))
        if ierror:
            return None
        # If imp.find_module returned an open file to rank 0, then we should
        # open the corresponding file for this process too.
        if stuff and stuff[1]:
            fp = open(pathname,stuff[1])

    try:
        m = imp.load_module(fqname, fp, pathname, stuff)
    finally:
        if fp: fp.close()
    if parent:
        setattr(parent, partname, m)
    return m


# The remaining functions are taken unmodified (except for the names)
# from knee.py.
def __determine_parent__(globals, level):
    if not globals or "__name__" not in globals:
        return None
    pname = globals['__name__']
    if "__path__" in globals:
        parent = sys.modules[pname]
        assert globals is parent.__dict__
        return parent
    if '.' in pname:
        if level > 0:
            end = len(pname)
            for l in range(level):
                i = pname.rfind('.', 0, end)
                end = i
        else:
            i = pname.rfind('.')
        pname = pname[:i]
        parent = sys.modules[pname]
        assert parent.__name__ == pname
        return parent
    return None

def __find_head_package__(parent, name):
    if '.' in name:
        i = name.find('.')
        head = name[:i]
        tail = name[i+1:]
    else:
        head = name
        tail = ""
    if parent:
        qname = "%s.%s" % (parent.__name__, head)
    else:
        qname = head
    q = __import_module__(head, qname, parent)
    if q: return q, tail
    if parent:
        qname = head
        parent = None
        q = __import_module__(head, qname, parent)
        if q: return q, tail
    raise ImportError("No module named " + qname)

def __load_tail__(q, tail):
    m = q
    while tail:
        i = tail.find('.')
        if i < 0: i = len(tail)
        head, tail = tail[:i], tail[i+1:]
        mname = "%s.%s" % (m.__name__, head)
        m = __import_module__(head, mname, m)
        if not m:
            raise ImportError("No module named " + mname)
    return m

def __ensure_fromlist__(m, fromlist, recursive=0):
    for sub in fromlist:
        if sub == "*":
            if not recursive:
                try:
                    all = m.__all__
                except AttributeError:
                    pass
                else:
                    __ensure_fromlist__(m, all, 1)
            continue
        if sub != "*" and not hasattr(m, sub):
            subname = "%s.%s" % (m.__name__, sub)
            submod = __import_module__(sub, subname, m)
            if not submod:
                raise ImportError("No module named " + subname)

# Now we import all the yt.mods items.
with mpi_import():
    if MPI.COMM_WORLD.rank == 0: print("Beginning parallel import block.")
    from yt.mods import *  # NOQA
    if MPI.COMM_WORLD.rank == 0: print("Ending parallel import block.")
