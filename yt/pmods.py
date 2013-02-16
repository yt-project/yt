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


"""This is an initial implementation of the finder/loader discussed at:
http://mail.scipy.org/pipermail/numpy-discussion/2012-March/061160.html

This is intended to take the place of MPI_Import.py. This version has
only been tested minimally, and is being made available primarily for
testing and preliminary benchmarking.

Known issues:
- Modules loaded via the Windows registry may be incorrectly hidden by
  a module of the same name in sys.path.
- If a file is added to a directory on sys.path, it won't be cached, so
  there may be precedence issues. If a file disappears or its permissions
  change, the import will fail.

Update (3/16/12): I've merged in a new version, simple_finder, described
below.

To use the finder, start a script off with the following:

import sys
from cached_import import finder
sys.meta_path.append(finder())

There are also variants of the finder that use MPI. The rank 0 process
builds the cache and then broadcasts it. For these, replace finder
with either pympi_finder or mpi4py_finder.

This finder works by building a cache mapping module names to
locations. The expensive parts of this process are the calls that
result in a stat. For that reason, we don't, by default, check whether
a module file is readable.

Since calls like os.isfile are expensive, I've added an alternate
version called simple_finder. Instead of figuring out where all of the
modules in sys.path are located, we just cache the contents of
directories on sys.path and use the standard probing algorithm for the
imports. This is much cheaper at startup and easier to maintain. It
appears to be a bit faster than the MPI-enabled finders, though that
will depend on the number of modules in sys.path as well as the number
of modules actually imported.
"""

import sys,os,imp

class finder(object):
    def __init__(self,skip_checks=True,build=True):
        """Build a finder object.

        Arguments:
        - skip_checks: Don't test whether modules are readable while building
                       the cache. This improves performace, but can cause an
                       unreadable file that looks like a Python module to
                       shadow a readable module with the same name later
                       in sys.path.
        -build: if set, build the cache now. This is used in the mpi4py_finder
                and pympi_finder extensions
        """
        # Store some suffix and module description information
        t = imp.get_suffixes()
        self.skip_checks = skip_checks
        self._suffixes = [x[0] for x in t] # in order of precedence
        self._rsuffixes = self._suffixes[::-1] # and in reverse order
        self._suffix_tuples = dict((x[0],tuple(x)) for x in t)

        # We store the value of sys.path in _syspath so we can keep track
        # of changes. _cache is a dictionary mapping module names to tuples
        # containing the information needed to load the module (path and
        # module description).
        if build:
            self._syspath = list(sys.path)
            self._build_cache()
        else: # For some subclasses
            self._syspath = []
            self._cache = {}

    def _build_cache(self):
        """Traverse sys.path, building (or re-building) the cache."""
        import os
        self._cache = {}
        for d in self._syspath:
            self._process_dir(os.path.realpath(d))

    def find_module(self,fullname,path=None):
        """Return self if 'fullname' is in sys.path (and isn't a builtin or
        frozen module)."""
        # First, make sure our cache is up-to-date. (We could combine
        # the append/prepend cases and more generally handle the case where
        # self._syspath is a sublist of the new sys.path, but is that worth
        # the effort? It's only beneficial if we encounter code where sys.path
        # is both prepended to and appended to, and there isn't an import
        # statement in between.
        if sys.path != self._syspath:
            stored_length = len(self._syspath)
            real_length = len(sys.path)
            rebuild = False
            # If sys.path isn't bigger, we need to rebuild the cache
            # but not before we update self._syspath.
            if real_length <= stored_length:
                rebuild = True
            # Some directories were prepended to the path, so add them.
            elif self._syspath == sys.path[-stored_length:]:
                for d in sys.path[real_length-stored_length-1::-1]:
                    self._process_dir(os.path.realpath(d),prepend=True)
            # Directories appended to the path.
            elif self._syspath == sys.path[:len(self._syspath)]:
                for d in sys.path[stored_length-real_length:]:
                    self._process_dir(os.path.realpath(d))
            # Path otherwise modified, so we need to rebuild the cache.
            else:
                rebuild = True

            # Now update self._syspath
            self._syspath = list(sys.path)
            if rebuild:
                self._build_cache()
            
        # Don't override builtin/frozen modules. TODO: Windows registry?
        if (fullname not in sys.builtin_module_names and
            not imp.is_frozen(fullname) and
            fullname in self._cache):
            #print "__IMPORTING ",fullname
            return self
        return None

    def load_module(self,fullname):
        """Load the module fullname using cached path."""
        if fullname in self._cache:
            if fullname in sys.modules:
                return sys.modules[fullname]
            pathname,desc = self._cache[fullname]
            #print "__LOADING ",fullname,pathname
            if os.path.isfile(pathname):
                # (If we're loading a PY_SOURCE file, the interpreter will
                # automatically check for a compiled (.py[c|o]) file.)
                with open(pathname,desc[1]) as f:
                    mod = imp.load_module(fullname,f,pathname,desc)
            # Not a file, so it's a package directory
            else:
                mod = imp.load_module(fullname,None,pathname,desc)
            mod.__loader__ = self
            return mod
        raise ImportError("This shouldn't happen!")


    # Build up a dict of modules (including package directories) found in a
    # directory. If this directory has been prepended to the path, we need to
    # overwrite any conflicting entries in the cache. To make sure precedence
    # is correct, we'll reverse the list of suffixes when we're prepending.
    #
    # Rather than add a lot of checks here to make sure we don't stomp on a
    # builtin module, we'll just reject these in find_module
    def _process_dir(self,dir,parent=None,prepend=False,visited=None):
        """Process a directory dir, looking for valid modules.

        Arguments:
        dir -- (an absolute, real path to a directory)
        parent -- parent module, in the case where dir is a package directory
        prepend -- True if dir has just been prepended to sys.path. In that
                   case, we'll replace existing cached entries with the same
                   module name.
        visited -- list of the real paths of visited directories. Used to
                   prevent infinite recursion in the case of symlink cycles
                   in package subdirectories.
        """
        import stat
        
        # Avoid symlink cycles in a package.
        if not visited:
            visited = [dir]
        elif dir not in visited:
            visited.append(dir)
        else:
            return

        # All files and subdirs. Store the name and the path.
        try:
            contents = dict((x,os.path.join(dir,x))
                            for x in os.listdir(dir))
        # Unreadable directory, so skip
        except OSError:
            return

        # If this is a possible package directory with no __init__.py, bail
        # out. If __init__.py is there, we need to see if there's an exising
        # module by that name. 
        if parent:
            if "__init__.py" not in contents:
                return
            if not (self.skip_checks or
                    os.access(os.path.join(dir,"__init__.py"),os.R_OK)):
                return
            if parent in self._cache and not prepend:
                return
            # Okay, this is a valid, non-duplicate module.
            self._cache[parent] = (dir,('','',imp.PKG_DIRECTORY))
            
        # Split contents into files & subdirs (only stat each one once)
        files = {}
        subdirs = {}
        for entry in contents:
            try:
                mode = os.stat(contents[entry]).st_mode
            except OSError:
                continue # couldn't read!
            if stat.S_ISDIR(mode) and (self.skip_checks or
                                       os.access(contents[entry],os.R_OK)):
                subdirs[entry] = contents[entry]
            elif stat.S_ISREG(mode) and (self.skip_checks or
                                         os.access(contents[entry],os.R_OK)):
                files[entry] = contents[entry]

        # Package directories have the highest precedence. But when prepend is
        # True, we need to reverse the order here. We'll do this with these
        # nested functions.
        def process_subdirs():
            for d in subdirs:
                fqname = parent+"."+d if parent else d # fully qualified name
                self._process_dir(os.path.join(dir,d),fqname,prepend,visited)

        def process_files():
            ordered_suffixes = self._rsuffixes if prepend else self._suffixes
            for s in ordered_suffixes:
                l = len(s)
                for f in files:
                    # Check for matching suffix.
                    if f[-l:] == s:
                        fqname = parent+"."+f[:-l] if parent else f[:-l]
                        if fqname not in self._cache or prepend:
                                self._cache[fqname] = (files[f],
                                                       self._suffix_tuples[s])

        if prepend:
            process_files()
            process_subdirs()
        else:
            process_subdirs()
            process_files()

                                
"""Finder that lets one MPI process do all of the initial caching.
"""
class pympi_finder(finder):        
    def __init__(self,skip_checks=True):
        import mpi
        if mpi.rank == 0:
            finder.__init__(self,skip_checks)
        else:
            finder.__init__(self,skip_checks,False)
        self._syspath,self._cache = mpi.bcast((self._syspath,self._cache))

"""Finder that lets one MPI process do all of the initial caching.
"""
class mpi4py_finder(finder):        
    def __init__(self,skip_checks=True):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            finder.__init__(self,skip_checks)
        else:
            finder.__init__(self,skip_checks,False)
        self._syspath,self._cache = comm.bcast((self._syspath,self._cache))

"""
Alternate version of cached_import. Instead of caching locations,
just cache directory contents. Then mimic the standard import probing
algorithm.

This has not been thoroughly tested!
"""

class simple_finder(object):    
    def __init__(self):
        # _contents is a representation of the files located in
        # sys.path (including, in the case of module packages, any
        # subdirectories that are encountered in the import process).
        # For each string in sys.path or subdirectory visited,
        # _contents contains a dict mapping the filenames in the
        # directory to full paths. If the string doesn't represent a
        # valid directory, then the dict is empty.
        self._contents = {}
        for d in sys.path:
            self._process_dir(d)

    # Search for a module 'name' in the cached directory listing for 'path'
    def _search(self,name,path):
        # If we haven't cached the directory, do so now.
        if path not in self._contents:
            self._process_dir(path)
        listing = self._contents[path]
        # First check for a package directory.
        try:
            if (name in listing and
                os.path.isdir(listing[name]) and
                os.path.isfile(os.path.join(listing[name],
                                            "__init__.py"))):
                return listing[name],('','',imp.PKG_DIRECTORY)
        except OSError:
            pass
        # Now check probe for name.so, namemodule.so, name.py, etc.
        for suffix in imp.get_suffixes():
            s = name+suffix[0]
            if (s in listing and
                os.path.isfile(listing[s]) and
                os.access(listing[s],os.R_OK)):
                return listing[s],suffix
        return None,None
    
    # Don't use this directly. We need more state than the load_module
    # signature allows, so we'll return a loader object for any module
    # that we have found.
    class _loader(object):
        def __init__(self,fullname,path,desc,finder):
            self._fullname = fullname
            self._path = path
            self._desc = desc
            self._finder = finder

        def load_module(self,fullname):
            """Load the module fullname using cached path."""
            if fullname != self._fullname:
                raise ImportError 
            if os.path.isfile(self._path): # check desc instead?
                with open(self._path,self._desc[1]) as f:
                    mod = imp.load_module(fullname,f,self._path,self._desc)
            # Not a file, so it's a package directory
            else:
                mod = imp.load_module(fullname,None,self._path,self._desc)
            mod.__loader__ = self._finder
            return mod

    # "Loader" for modules that have already been imported.
    class _null_loader(object):
        def load_module(self,fullname):
            return sys.modules[fullname]

    def find_module(self,fullname,path=None):
        """Return self if 'fullname' is in sys.path (and isn't a builtin or
        frozen module)."""
        if fullname in sys.modules:
            return simple_finder._null_loader()
        # Don't override builtin/frozen modules. TODO: Windows registry?
        if (fullname not in sys.builtin_module_names and
            not imp.is_frozen(fullname)):
            if path:
                iterpath = path
                name = fullname.split('.')[-1]
            else:
                iterpath = sys.path
                name = fullname
            for dir in iterpath:
                loadpath,desc = self._search(name,dir)
                if loadpath:
                    break
            #print "__IMPORTING ",fullname
            if loadpath:
                return simple_finder._loader(fullname,loadpath,desc,self)
        return None

    def _process_dir(self,dir):
        """
        Arguments:
        dir -- 
        """
        # All files and subdirs. Store the name and the path.
        try:
            contents = dict((x,os.path.join(dir,x))
                            for x in os.listdir(dir))
        # Unreadable directory, so skip
        except OSError:
            contents = {}

        self._contents[dir] = contents

# Now we import all the yt.mods items.
import sys
sys.meta_path.append(mpi4py_finder())
from yt.mods import *
