#!/usr/bin/env python
import setuptools
import os, sys, os.path, glob

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

def check_for_freetype():
    # First up: environment
    if "FTYPE_DIR" in os.environ:
        freetype_dir = os.environ["FTYPE_DIR"]
        freetype_inc = os.path.join(freetype_dir, "include")
        freetype_lib = os.path.join(freetype_dir, "lib")
        print "FTYPE_LOCATION: FTYPE_DIR: %s, %s" % (freetype_inc, freetype_lib)
        return (freetype_inc, freetype_lib)
    # Next up, we try freetype.cfg
    elif os.path.exists("freetype.cfg"):
        freetype_dir = open("freetype.cfg").read().strip()
        freetype_inc = os.path.join(freetype_dir, "include")
        freetype_lib = os.path.join(freetype_dir, "lib")
        print "FTYPE_LOCATION: freetype.cfg: %s, %s" % (freetype_inc, freetype_lib)
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
               os.path.isfile(os.path.join(freetype_dir, "include", "ft2build.h")):
                freetype_inc = os.path.join(freetype_dir, "include")
                freetype_lib = os.path.join(freetype_dir, "lib")
                print "FTYPE_LOCATION: freetype found in: %s, %s" % (freetype_inc, freetype_lib)
                return freetype_inc, freetype_lib
    except ImportError:
        pass
    # X11 is where it's located by default on OSX, although I am slightly
    # reluctant to link against that one.
    for freetype_dir in ["/usr/", "/usr/local/", "/usr/X11/"]:
        if os.path.isfile(os.path.join(freetype_dir, "include", "ft2build.h")):
            if os.path.isdir(os.path.join(freetype_dir, "include")) and \
               os.path.isfile(os.path.join(freetype_dir, "include", "ft2build.h")):
                freetype_inc = os.path.join(freetype_dir, "include")
                freetype_lib = os.path.join(freetype_dir, "lib")
                print "FTYPE_LOCATION: freetype found in: %s, %s" % (freetype_inc, freetype_lib)
                return freetype_inc, freetype_lib
    print "Reading freetype location from freetype.cfg failed."
    print "Please place the base directory of your freetype install in freetype.cfg and restart."
    print "(ex: \"echo '/usr/local/' > freetype.cfg\" )"
    print "You can locate this by looking for the file ft2build.h"
    sys.exit(1)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lib',parent_package,top_path)
    png_inc, png_lib = check_for_png()
    freetype_inc, freetype_lib = check_for_freetype()
    # Because setjmp.h is included by lots of things, and because libpng hasn't
    # always properly checked its header files (see
    # https://bugzilla.redhat.com/show_bug.cgi?id=494579 ) we simply disable
    # support for setjmp.
    config.add_extension("CICDeposit", 
                ["yt/utilities/lib/CICDeposit.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("ContourFinding", 
                ["yt/utilities/lib/ContourFinding.pyx",
                 "yt/utilities/lib/union_find.c"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("DepthFirstOctree", 
                ["yt/utilities/lib/DepthFirstOctree.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("fortran_reader", 
                ["yt/utilities/lib/fortran_reader.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("freetype_writer", 
                ["yt/utilities/lib/freetype_writer.pyx"],
                include_dirs = [freetype_inc,os.path.join(freetype_inc, "freetype2"),
                                "yt/utilities/lib"],
                library_dirs = [freetype_lib], libraries=["freetype"],
                depends=["yt/utilities/lib/freetype_includes.h"])
    config.add_extension("geometry_utils", 
                ["yt/utilities/lib/geometry_utils.pyx"],
               extra_compile_args=['-fopenmp'],
               extra_link_args=['-fopenmp'],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("Interpolators", 
                ["yt/utilities/lib/Interpolators.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("marching_cubes", 
                ["yt/utilities/lib/marching_cubes.pyx",
                 "yt/utilities/lib/FixedInterpolator.c"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/utilities/lib/fixed_interpolator.pxd",
                         "yt/utilities/lib/FixedInterpolator.h",
                ])
    config.add_extension("misc_utilities", 
                ["yt/utilities/lib/misc_utilities.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("Octree", 
                ["yt/utilities/lib/Octree.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("png_writer", 
                ["yt/utilities/lib/png_writer.pyx"],
                define_macros=[("PNG_SETJMP_NOT_SUPPORTED", True)],
                include_dirs=[png_inc],
                library_dirs=[png_lib],
                libraries=["m", "png"],
                depends=["yt/utilities/lib/fp_utils.pxd"]),
    config.add_extension("PointsInVolume", 
                ["yt/utilities/lib/PointsInVolume.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("QuadTree", 
                ["yt/utilities/lib/QuadTree.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("RayIntegrators", 
                ["yt/utilities/lib/RayIntegrators.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("VolumeIntegrator", 
               ["yt/utilities/lib/VolumeIntegrator.pyx",
                "yt/utilities/lib/FixedInterpolator.c",
                "yt/utilities/lib/kdtree.c"] +
                 glob.glob("yt/utilities/lib/healpix_*.c"), 
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               depends = ["yt/utilities/lib/VolumeIntegrator.pyx",
                          "yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/healpix_interface.pxd",
                          "yt/utilities/lib/endian_swap.h",
                          "yt/utilities/lib/FixedInterpolator.h",
                          "yt/utilities/lib/healpix_vectors.h",
                          "yt/utilities/lib/kdtree.h",
                          "yt/utilities/lib/healpix_ang2pix_nest.c",
                          "yt/utilities/lib/healpix_mk_pix2xy.c",
                          "yt/utilities/lib/healpix_mk_xy2pix.c",
                          "yt/utilities/lib/healpix_pix2ang_nest.c",
                          "yt/utilities/lib/healpix_pix2vec_nest.c",
                          "yt/utilities/lib/healpix_vec2pix_nest.c"]
          )
    config.add_extension("grid_traversal", 
               ["yt/utilities/lib/grid_traversal.pyx",
                "yt/utilities/lib/FixedInterpolator.c",
                "yt/utilities/lib/kdtree.c"] +
                 glob.glob("yt/utilities/lib/healpix_*.c"), 
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               extra_compile_args=['-fopenmp'],
               extra_link_args=['-fopenmp'],
               depends = ["yt/utilities/lib/VolumeIntegrator.pyx",
                          "yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/kdtree.h",
                          "yt/utilities/lib/FixedInterpolator.h",
                          "yt/utilities/lib/fixed_interpolator.pxd",
                          "yt/utilities/lib/field_interpolation_tables.pxd",
                          ]
          )
    if os.environ.get("GPERFTOOLS", "no").upper() != "NO":
        gpd = os.environ["GPERFTOOLS"]
        idir = os.path.join(gpd, "include")
        ldir = os.path.join(gpd, "lib")
        print "INCLUDE AND LIB DIRS", idir, ldir
        config.add_extension("perftools_wrap",
                ["yt/utilities/lib/perftools_wrap.pyx"],
                libraries=["profiler"],
                library_dirs = [ldir],
                include_dirs = [idir],
            )
    config.make_config_py() # installs __config__.py
    return config
