#!/usr/bin/env python
import setuptools
import os, sys

def check_for_hdf5():
    # First up: HDF5_DIR in environment
    if "HDF5_DIR" in os.environ:
        hdf5_dir = os.environ["HDF5_DIR"]
        hdf5_inc = os.path.join(hdf5_dir, "include")
        hdf5_lib = os.path.join(hdf5_dir, "lib")
        print "HDF5_LOCATION: HDF5_DIR: %s, %s" % (hdf5_inc, hdf5_lib)
        return (hdf5_inc, hdf5_lib)
    # Next up, we try hdf5.cfg
    elif os.path.exists("hdf5.cfg"):
        hdf5_dir = open("hdf5.cfg").read().strip()
        hdf5_inc = os.path.join(hdf5_dir, "include")
        hdf5_lib = os.path.join(hdf5_dir, "lib")
        print "HDF5_LOCATION: hdf5.cfg: %s, %s" % (hdf5_inc, hdf5_lib)
        return (hdf5_inc, hdf5_lib)
    # Now we see if ctypes can help us:
    try:
        import ctypes.util
        hdf5_libfile = ctypes.util.find_library("hdf5")
        if hdf5_libfile is not None and os.path.isfile(hdf5_libfile):
            # Now we've gotten a library, but we'll need to figure out the
            # includes if this is going to work.  It feels like there is a
            # better way to pull off two directory names.
            hdf5_dir = os.path.dirname(os.path.dirname(hdf5_libfile))
            if os.path.isdir(os.path.join(hdf5_dir, "include")) and \
               os.path.isfile(os.path.join(hdf5_dir, "include", "hdf5.h")):
                hdf5_inc = os.path.join(hdf5_dir, "include")
                hdf5_lib = os.path.join(hdf5_dir, "lib")
                print "HDF5_LOCATION: HDF5 found in: %s, %s" % (hdf5_inc, hdf5_lib)
                return hdf5_inc, hdf5_lib
    except ImportError:
        pass
    print "Reading HDF5 location from hdf5.cfg failed."
    print "Please place the base directory of your HDF5 install in hdf5.cfg and restart."
    print "(ex: \"echo '/usr/local/' > hdf5.cfg\" )"
    sys.exit(1)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lagos',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.make_svn_version_py()
    config.add_extension("PointCombine", "yt/lagos/PointCombine.c", libraries=["m"])
    config.add_subpackage("hop")
    config.add_subpackage("fof")
    config.add_subpackage("parallelHOP")
    hdf5_inc, hdf5_lib = check_for_hdf5()
    include_dirs=[hdf5_inc]
    library_dirs=[hdf5_lib]
    config.add_extension("HDF5LightReader", "yt/lagos/HDF5LightReader.c",
                         define_macros=[("H5_USE_16_API",True)],
                         libraries=["m","hdf5"],
                         library_dirs=library_dirs, include_dirs=include_dirs)
    return config
