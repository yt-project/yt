#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob


# snatched from PyTables
def add_from_path(envname, dirs):
    try:
        dirs.extend(os.environ[envname].split(os.pathsep))
    except KeyError:
        pass


# snatched from PyTables
def add_from_flags(envname, flag_key, dirs):
    for flag in os.environ.get(envname, "").split():
        if flag.startswith(flag_key):
            dirs.append(flag[len(flag_key):])


# snatched from PyTables
def get_default_dirs():
    default_header_dirs = []
    add_from_path("CPATH", default_header_dirs)
    add_from_path("C_INCLUDE_PATH", default_header_dirs)
    add_from_flags("CPPFLAGS", "-I", default_header_dirs)
    default_header_dirs.extend(['/usr/include', '/usr/local/include'])

    default_library_dirs = []
    add_from_flags("LDFLAGS", "-L", default_library_dirs)
    default_library_dirs.extend(
        os.path.join(_tree, _arch)
        for _tree in ('/', '/usr', '/usr/local')
        for _arch in ('lib64', 'lib')
    )
    return default_header_dirs, default_library_dirs


def get_hdf5_include(header_dirs):
    for inc_prefix in header_dirs:
        if os.path.isfile(os.path.join(inc_prefix, "hdf5.h")):
            return inc_prefix
    return None


def get_hdf5_lib(lib_dirs):
    import ctypes
    hdf5_libfile = ctypes.util.find_library("hdf5")
    if os.path.isfile(hdf5_libfile):
        return os.path.dirname(hdf5_libfile)
    for lib_dir in lib_dirs:
        try:
            ctypes.CDLL(os.path.join(lib_dir, hdf5_libfile))
            return lib_dir
        except OSError:
            pass
    return None


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
    if os.name == 'posix':
        default_header_dirs, default_library_dirs = get_default_dirs()
        hdf5_inc = get_hdf5_include(default_header_dirs)
        hdf5_lib = get_hdf5_lib(default_library_dirs)
        if None not in (hdf5_inc, hdf5_lib):
            print(
                "HDF5_LOCATION: HDF5 found in: %s, %s" % (hdf5_inc, hdf5_lib))
            return (hdf5_inc, hdf5_lib)

    # Now we see if ctypes can help us on non posix platform
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
                print "HDF5_LOCATION: HDF5 found in: %s, %s" % (hdf5_inc,
                                                                hdf5_lib)
                return (hdf5_inc, hdf5_lib)
    except ImportError:
        pass
    print "Reading HDF5 location from hdf5.cfg failed."
    print "Please place the base directory of your"
    print "HDF5 install in hdf5.cfg and restart."
    print "(ex: \"echo '/usr/local/' > hdf5.cfg\" )"
    sys.exit(1)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('utilities', parent_package, top_path)
    config.add_subpackage("amr_kdtree")
    config.add_subpackage("poster")
    config.add_subpackage("answer_testing")
    config.add_subpackage("delaunay")  # From SciPy, written by Robert Kern
    config.add_subpackage("kdtree")
    config.add_subpackage("spatial")
    config.add_subpackage("grid_data_format")
    config.add_subpackage("parallel_tools")
    config.add_subpackage("lib")
    config.add_extension("data_point_utilities",
                         "yt/utilities/data_point_utilities.c",
                         libraries=["m"])
    config.add_subpackage("tests")
    hdf5_inc, hdf5_lib = check_for_hdf5()
    include_dirs = [hdf5_inc]
    library_dirs = [hdf5_lib]
    config.add_extension("hdf5_light_reader",
                         "yt/utilities/hdf5_light_reader.c",
                         define_macros=[("H5_USE_16_API", True)],
                         libraries=["m", "hdf5"],
                         library_dirs=library_dirs, include_dirs=include_dirs)
    config.add_extension("libconfig_wrapper",
                         ["yt/utilities/libconfig_wrapper.pyx"] +
                         glob.glob("yt/utilities/_libconfig/*.c"),
                         include_dirs=["yt/utilities/_libconfig/"],
                         define_macros=[("HAVE_XLOCALE_H", True)]
                         )
    config.make_config_py()  # installs __config__.py
    # config.make_svn_version_py()
    return config
