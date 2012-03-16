#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob


def check_for_png():
    # First up: HDF5_DIR in environment
    if "PNG_DIR" in os.environ:
        png_dir = os.environ["PNG_DIR"]
        png_inc = os.path.join(png_dir, "include")
        png_lib = os.path.join(png_dir, "lib")
        print "PNG_LOCATION: PNG_DIR: %s, %s" % (png_inc, png_lib)
        return (png_inc, png_lib)
    # Next up, we try png.cfg
    elif os.path.exists("png.cfg"):
        png_dir = open("png.cfg").read().strip()
        png_inc = os.path.join(png_dir, "include")
        png_lib = os.path.join(png_dir, "lib")
        print "PNG_LOCATION: png.cfg: %s, %s" % (png_inc, png_lib)
        return (png_inc, png_lib)
    # Now we see if ctypes can help us:
    try:
        import ctypes.util
        png_libfile = ctypes.util.find_library("png")
        if png_libfile is not None and os.path.isfile(png_libfile):
            # Now we've gotten a library, but we'll need to figure out the
            # includes if this is going to work.  It feels like there is a
            # better way to pull off two directory names.
            png_dir = os.path.dirname(os.path.dirname(png_libfile))
            if os.path.isdir(os.path.join(png_dir, "include")) and \
               os.path.isfile(os.path.join(png_dir, "include", "png.h")):
                png_inc = os.path.join(png_dir, "include")
                png_lib = os.path.join(png_dir, "lib")
                print "PNG_LOCATION: png found in: %s, %s" % (png_inc, png_lib)
                return png_inc, png_lib
    except ImportError:
        pass
    # X11 is where it's located by default on OSX, although I am slightly
    # reluctant to link against that one.
    for png_dir in ["/usr/", "/usr/local/", "/usr/X11/"]:
        if os.path.isfile(os.path.join(png_dir, "include", "png.h")):
            if os.path.isdir(os.path.join(png_dir, "include")) and \
               os.path.isfile(os.path.join(png_dir, "include", "png.h")):
                png_inc = os.path.join(png_dir, "include")
                png_lib = os.path.join(png_dir, "lib")
                print "PNG_LOCATION: png found in: %s, %s" % (png_inc, png_lib)
                return png_inc, png_lib
    print "Reading png location from png.cfg failed."
    print "Please place the base directory of your png"
    print "install in png.cfg and restart."
    print "(ex: \"echo '/usr/local/' > png.cfg\" )"
    sys.exit(1)


def check_for_freetype():
    # First up: environment
    if "FTYPE_DIR" in os.environ:
        freetype_dir = os.environ["FTYPE_DIR"]
        freetype_inc = os.path.join(freetype_dir, "include")
        freetype_lib = os.path.join(freetype_dir, "lib")
        print "FTYPE_LOCATION: FTYPE_DIR: %s, %s" % (freetype_inc,
            freetype_lib)
        return (freetype_inc, freetype_lib)
    # Next up, we try freetype.cfg
    elif os.path.exists("freetype.cfg"):
        freetype_dir = open("freetype.cfg").read().strip()
        freetype_inc = os.path.join(freetype_dir, "include")
        freetype_lib = os.path.join(freetype_dir, "lib")
        print "FTYPE_LOCATION: freetype.cfg: %s, %s" % (freetype_inc,
            freetype_lib)
        return (freetype_inc, freetype_lib)
    # Now we see if ctypes can help us:
    try:
        import ctypes.util
        freetype_libfile = ctypes.util.find_library("freetype")
        if freetype_libfile is not None and os.path.isfile(freetype_libfile):
            # Now we've gotten a library, but we'll need to figure out the
            # includes if this is going to work.  It feels like there is a
            # better way to pull off two directory names.
            freetype_dir = os.path.dirname(os.path.dirname(freetype_libfile))
            if os.path.isdir(os.path.join(freetype_dir, "include")) and \
                os.path.isfile(os.path.join(freetype_dir, "include",
                    "ft2build.h")):
                freetype_inc = os.path.join(freetype_dir, "include")
                freetype_lib = os.path.join(freetype_dir, "lib")
                print "FTYPE_LOCATION: freetype found in: %s, %s" % \
                    (freetype_inc, freetype_lib)
                return freetype_inc, freetype_lib
    except ImportError:
        pass
    # X11 is where it's located by default on OSX, although I am slightly
    # reluctant to link against that one.
    for freetype_dir in ["/usr/", "/usr/local/", "/usr/X11/"]:
        if os.path.isfile(os.path.join(freetype_dir, "include", "ft2build.h")):
            if os.path.isdir(os.path.join(freetype_dir, "include")) and \
                os.path.isfile(os.path.join(freetype_dir, "include",
                    "ft2build.h")):
                freetype_inc = os.path.join(freetype_dir, "include")
                freetype_lib = os.path.join(freetype_dir, "lib")
                print "FTYPE_LOCATION: freetype found in: %s, %s" % \
                    (freetype_inc, freetype_lib)
                return freetype_inc, freetype_lib
    print "Reading freetype location from freetype.cfg failed."
    print "Please place the base directory of your freetype"
    print "install in freetype.cfg and restart."
    print "(ex: \"echo '/usr/local/' > freetype.cfg\" )"
    print "You can locate this by looking for the file ft2build.h"
    sys.exit(1)


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
                print "HDF5_LOCATION: HDF5 found in: %s, %s" % (hdf5_inc,
                    hdf5_lib)
                return hdf5_inc, hdf5_lib
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
    png_inc, png_lib = check_for_png()
    freetype_inc, freetype_lib = check_for_freetype()
    # Because setjmp.h is included by lots of things,
    # and because libpng hasn't
    # always properly checked its header files (see
    # https://bugzilla.redhat.com/show_bug.cgi?id=494579 ) we simply disable
    # support for setjmp.
    config.add_subpackage("amr_kdtree")
    config.add_subpackage("answer_testing")
    config.add_subpackage("delaunay")  # From SciPy, written by Robert Kern
    config.add_subpackage("kdtree")
    config.add_data_files(('kdtree', ['kdtree/fKDpy.so']))
    config.add_subpackage("spatial")
    config.add_subpackage("parallel_tools")
    config.add_extension("data_point_utilities",
                "yt/utilities/data_point_utilities.c", libraries=["m"])
    hdf5_inc, hdf5_lib = check_for_hdf5()
    include_dirs = [hdf5_inc]
    library_dirs = [hdf5_lib]
    config.add_extension("hdf5_light_reader",
                        "yt/utilities/hdf5_light_reader.c",
                         define_macros=[("H5_USE_16_API", True)],
                         libraries=["m", "hdf5"],
                         library_dirs=library_dirs, include_dirs=include_dirs)
    config.add_extension("amr_utils",
        ["yt/utilities/amr_utils.pyx",
         "yt/utilities/_amr_utils/FixedInterpolator.c",
         "yt/utilities/_amr_utils/kdtree.c",
         "yt/utilities/_amr_utils/union_find.c"] +
         glob.glob("yt/utilities/_amr_utils/healpix_*.c"),
        define_macros=[("PNG_SETJMP_NOT_SUPPORTED", True)],
        include_dirs=["yt/utilities/_amr_utils/", png_inc,
                      freetype_inc, os.path.join(freetype_inc, "freetype2")],
        library_dirs=[png_lib, freetype_lib],
        libraries=["m", "png", "freetype"],
        depends=glob.glob("yt/utilities/_amr_utils/*.pyx") +
                glob.glob("yt/utilities/_amr_utils/*.h") +
                glob.glob("yt/utilities/_amr_utils/*.c"),
        )
    #config.add_extension("voropp",
    #    ["yt/utilities/voropp.pyx"],
    #    language="c++",
    #    include_dirs=["yt/utilities/voro++"])
    config.add_extension("libconfig_wrapper",
        ["yt/utilities/libconfig_wrapper.pyx"] +
         glob.glob("yt/utilities/_libconfig/*.c"),
        include_dirs=["yt/utilities/_libconfig/"],
        define_macros=[("HAVE_XLOCALE_H", True)]
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
