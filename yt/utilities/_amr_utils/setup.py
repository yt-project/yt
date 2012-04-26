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
                ["yt/utilities/_amr_utils/CICDeposit.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("ContourFinding", 
                ["yt/utilities/_amr_utils/ContourFinding.pyx",
                 "yt/utilities/_amr_utils/union_find.c"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("DepthFirstOctree", 
                ["yt/utilities/_amr_utils/DepthFirstOctree.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("fortran_reader", 
                ["yt/utilities/_amr_utils/fortran_reader.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("freetype_writer", 
                ["yt/utilities/_amr_utils/freetype_writer.pyx"],
                include_dirs = [freetype_inc,os.path.join(freetype_inc, "freetype2")],
                library_dirs = [freetype_lib], libraries=["freetype"],
                depends=["yt/utilities/_amr_utils/freetype_includes.h"])
    config.add_extension("geometry_utils", 
                ["yt/utilities/_amr_utils/geometry_utils.pyx"],
               extra_compile_args=['-fopenmp'],
               extra_link_args=['-fopenmp'],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("Interpolators", 
                ["yt/utilities/_amr_utils/Interpolators.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("marching_cubes", 
                ["yt/utilities/_amr_utils/marching_cubes.pyx",
                 "yt/utilities/_amr_utils/FixedInterpolator.c"],
                libraries=["m"],
                depends=["yt/utilities/_amr_utils/fp_utils.pxd",
                         "yt/utilities/_amr_utils/fixed_interpolator.pxd",
                         "yt/utilities/_amr_utils/FixedInterpolator.h",
                ])
    config.add_extension("misc_utilities", 
                ["yt/utilities/_amr_utils/misc_utilities.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("Octree", 
                ["yt/utilities/_amr_utils/Octree.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("png_writer", 
                ["yt/utilities/_amr_utils/png_writer.pyx"],
                define_macros=[("PNG_SETJMP_NOT_SUPPORTED", True)],
                include_dirs=[png_inc],
                library_dirs=[png_lib],
                libraries=["m", "png"],
                depends=["yt/utilities/_amr_utils/fp_utils.pxd"]),
    config.add_extension("PointsInVolume", 
                ["yt/utilities/_amr_utils/PointsInVolume.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("QuadTree", 
                ["yt/utilities/_amr_utils/QuadTree.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("RayIntegrators", 
                ["yt/utilities/_amr_utils/RayIntegrators.pyx"],
                libraries=["m"], depends=["yt/utilities/_amr_utils/fp_utils.pxd"])
    config.add_extension("VolumeIntegrator", 
               ["yt/utilities/_amr_utils/VolumeIntegrator.pyx",
                "yt/utilities/_amr_utils/FixedInterpolator.c",
                "yt/utilities/_amr_utils/kdtree.c"] +
                 glob.glob("yt/utilities/_amr_utils/healpix_*.c"), 
               include_dirs=["yt/utilities/_amr_utils/"],
               libraries=["m"], 
               depends = ["yt/utilities/_amr_utils/VolumeIntegrator.pyx",
                          "yt/utilities/_amr_utils/fp_utils.pxd",
                          "yt/utilities/_amr_utils/healpix_interface.pxd",
                          "yt/utilities/_amr_utils/endian_swap.h",
                          "yt/utilities/_amr_utils/FixedInterpolator.h",
                          "yt/utilities/_amr_utils/healpix_vectors.h",
                          "yt/utilities/_amr_utils/kdtree.h",
                          "yt/utilities/_amr_utils/healpix_ang2pix_nest.c",
                          "yt/utilities/_amr_utils/healpix_mk_pix2xy.c",
                          "yt/utilities/_amr_utils/healpix_mk_xy2pix.c",
                          "yt/utilities/_amr_utils/healpix_pix2ang_nest.c",
                          "yt/utilities/_amr_utils/healpix_pix2vec_nest.c",
                          "yt/utilities/_amr_utils/healpix_vec2pix_nest.c"]
          )
    config.add_extension("grid_traversal", 
               ["yt/utilities/_amr_utils/grid_traversal.pyx",
                "yt/utilities/_amr_utils/FixedInterpolator.c",
                "yt/utilities/_amr_utils/kdtree.c"] +
                 glob.glob("yt/utilities/_amr_utils/healpix_*.c"), 
               include_dirs=["yt/utilities/_amr_utils/"],
               libraries=["m"], 
               extra_compile_args=['-fopenmp'],
               extra_link_args=['-fopenmp'],
               depends = ["yt/utilities/_amr_utils/VolumeIntegrator.pyx",
                          "yt/utilities/_amr_utils/fp_utils.pxd",
                          "yt/utilities/_amr_utils/kdtree.h",
                          "yt/utilities/_amr_utils/FixedInterpolator.h",
                          "yt/utilities/_amr_utils/fixed_interpolator.pxd",
                          ]
          )
    if os.environ.get("GPERFTOOLS", "no").upper() != "NO":
        gpd = os.environ["GPERFTOOLS"]
        idir = os.path.join(gpd, "include")
        ldir = os.path.join(gpd, "lib")
        print "INCLUDE AND LIB DIRS", idir, ldir
        config.add_extension("perftools_wrap",
                ["yt/utilities/_amr_utils/perftools_wrap.pyx"],
                libraries=["profiler"],
                library_dirs = [ldir],
                include_dirs = [idir],
            )
    config.make_config_py() # installs __config__.py
    return config
