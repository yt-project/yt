#!/usr/bin/env python
import setuptools
import os, sys

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
    print "Please place the base directory of your png install in png.cfg and restart."
    print "(ex: \"echo '/usr/local/' > png.cfg\" )"
    sys.exit(1)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('yt', parent_package, top_path)
    config.add_subpackage('lagos')
    config.add_subpackage('raven')
    config.add_subpackage('fido')
    config.add_subpackage('reason')
    config.add_subpackage('extensions')
    config.add_subpackage('parallel_tools')
    config.add_subpackage('data_objects')
    png_inc, png_lib = check_for_png()
    include_dirs=[png_inc]
    library_dirs=[png_lib]
    # Because setjmp.h is included by lots of things, and because libpng hasn't
    # always properly checked its header files (see
    # https://bugzilla.redhat.com/show_bug.cgi?id=494579 ) we simply disable
    # support for setjmp.
    config.add_extension("amr_utils", 
        ["yt/amr_utils.c", "yt/_amr_utils/FixedInterpolator.c"], 
        define_macros=[("PNG_SETJMP_NOT_SUPPORTED", True)],
        include_dirs=["yt/_amr_utils/", png_inc],
        library_dirs=[png_lib],
        libraries=["m", "png"])
    config.make_config_py()
    config.make_svn_version_py()
    return config

if __name__ == '__main__':
    # Remove current working directory from sys.path
    # to avoid importing numpy.distutils as Python std. distutils:
    from numpy.distutils.core import setup
    setup(configuration=configuration)
