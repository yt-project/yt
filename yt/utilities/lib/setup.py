#!/usr/bin/env python
import setuptools
import os, sys, os.path, glob, \
    tempfile, subprocess, shutil, \
    platform

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
    default_library_dirs = []

    add_from_path("CPATH", default_header_dirs)
    add_from_path("C_INCLUDE_PATH", default_header_dirs)
    add_from_flags("CPPFLAGS", "-I", default_header_dirs)
    default_header_dirs.extend(
        ['/usr/include', '/usr/local/include', '/usr/X11']
    )

    _archs = ['lib64', 'lib']
    if platform.system() == 'Linux':
        distname, version, did = platform.linux_distribution()
        if distname.lower() in ('ubuntu', 'debian'):
            _archs.extend(
                ['lib/x86_64-linux-gnu',
                 'lib/i686-linux-gnu',
                 'lib/i386-linux-gnu']
            )

    add_from_flags("LDFLAGS", "-L", default_library_dirs)
    default_library_dirs.extend(
        os.path.join(_tree, _arch)
        for _tree in ('/usr', '/usr/local', '/usr/X11', '/')
        for _arch in _archs
    )
    return default_header_dirs, default_library_dirs


def get_location_from_env(env):
    env_dir = os.environ[env]
    env_inc = os.path.join(env_dir, "include")
    env_lib = os.path.join(env_dir, "lib")
    print("%s_LOCATION: %s: %s, %s"
          % (env.split('_')[0], env, env_inc, env_lib))
    return (env_inc, env_lib)


def get_location_from_cfg(cfg):
    cfg_dir = open(cfg).read().strip()
    cfg_inc = os.path.join(cfg_dir, "include")
    cfg_lib = os.path.join(cfg_dir, "lib")
    print("%s_LOCATION: %s: %s, %s"
          % (cfg.split('.')[0].upper(), cfg, cfg_inc, cfg_lib))
    return (cfg_inc, cfg_lib)


def check_prefix(inc_dir, lib_dir):
    if platform.system() == 'Linux':
        distname, version, did = platform.linux_distribution()
        if distname.lower() in ('ubuntu', 'debian'):
            print("Since you are using multiarch distro it's hard to detect")
            print("whether library matches the header file. We will assume")
            print("it does. If you encounter any build failures please use")
            print("proper cfg files to provide path to the dependencies")
            print("")
            return (inc_dir, lib_dir)
    prefix = os.path.commonprefix([inc_dir, lib_dir]).rstrip('/\\')
    if prefix is not '' and prefix == os.path.dirname(inc_dir):
        return (inc_dir, lib_dir)
    else:
        print("It seems that include prefix is different from lib prefix")
        print("Please use either env variable or cfg to set proper path")
        return (None, None)


def get_location_from_ctypes(header, library):
    yt_inst = os.environ.get('YT_DEST')
    if yt_inst is not None:
        # since we prefer installation via script, make sure
        # that YT_DEST path take precedence above all else
        return (os.path.join(yt_inst, 'include'), os.path.join(yt_inst, 'lib'))

    try:
        import ctypes
        import ctypes.util
    except ImportError:
        return (None, None)

    target_inc, target_libdir = None, None
    default_header_dirs, default_library_dirs = get_default_dirs()
    for inc_prefix in default_header_dirs:
        if os.path.isfile(os.path.join(inc_prefix, header)):
            target_inc = inc_prefix

    target_libfile = ctypes.util.find_library(library)
    if None in (target_inc, target_libfile):
        # either header or lib was not found, abort now
        return (None, None)
    if os.path.isfile(target_libfile):
        return check_prefix(target_inc, os.path.dirname(target_libfile))
    for lib_dir in default_library_dirs:
        try:
            ctypes.CDLL(os.path.join(lib_dir, target_libfile))
            target_libdir = lib_dir
        except OSError:
            pass
    return check_prefix(target_inc, target_libdir)


def check_for_dependencies(env, cfg, header, library):
    # First up: check in environment
    if env in os.environ:
        return get_location_from_env(env)
    # Next up, we try config file
    elif os.path.exists(cfg):
        return get_location_from_cfg(cfg)
    # Now we see if ctypes can help us
    if os.name == 'posix':
        target_inc, target_lib = get_location_from_ctypes(header, library)
    if None not in (target_inc, target_lib):
        print(
            "%s_LOCATION: %s found via ctypes in: %s, %s"
            % (env.split('_')[0], env.split('_')[0], target_inc, target_lib)
        )
        return (target_inc, target_lib)

    print("Reading %s location from %s failed." % (env.split('_')[0], cfg))
    print("Please place the base directory of your")
    print("%s install in %s and restart." % (env.split('_')[0], cfg))
    print("(ex: \"echo '/usr/local/' > %s\" )" % cfg)
    print("You can locate the path by looking for %s" % header)
    sys.exit(1)

def check_for_png():
    return check_for_dependencies("PNG_DIR", "png.cfg", "png.h", "png")


def check_for_openmp():
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')

        # Attempt to compile a test script.
        # See http://openmp.org/wp/openmp-compilers/
        filename = r'test.c'
        file = open(filename,'w', 0)
        file.write(
            "#include <omp.h>\n"
            "#include <stdio.h>\n"
            "int main() {\n"
            "#pragma omp parallel\n"
            "printf(\"Hello from thread %d, nthreads %d\\n\", omp_get_thread_num(), omp_get_num_threads());\n"
            "}"
            )
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call([compiler, '-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)

        # Clean up
        file.close()
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return exit_code == 0

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lib',parent_package,top_path)
    png_inc, png_lib = check_for_png()
    if check_for_openmp() == True:
        omp_args = ['-fopenmp']
    else:
        omp_args = None
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
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/utilities/lib/amr_kdtools.pxd"])
    config.add_extension("DepthFirstOctree", 
                ["yt/utilities/lib/DepthFirstOctree.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("fortran_reader", 
                ["yt/utilities/lib/fortran_reader.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("geometry_utils", 
                ["yt/utilities/lib/geometry_utils.pyx"],
               extra_compile_args=omp_args,
               extra_link_args=omp_args,
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("Interpolators", 
                ["yt/utilities/lib/Interpolators.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("alt_ray_tracers", 
                ["yt/utilities/lib/alt_ray_tracers.pyx"],
                libraries=["m"], depends=[])
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
                "yt/utilities/lib/kdtree.c"],
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               depends = ["yt/utilities/lib/VolumeIntegrator.pyx",
                          "yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/healpix_interface.pxd",
                          "yt/utilities/lib/endian_swap.h",
                          "yt/utilities/lib/FixedInterpolator.h",
                          "yt/utilities/lib/kdtree.h"],
          )
    config.add_extension("mesh_utilities",
              ["yt/utilities/lib/mesh_utilities.pyx"],
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               depends = ["yt/utilities/lib/fp_utils.pxd",
                          ],
          )
    config.add_extension("grid_traversal", 
               ["yt/utilities/lib/grid_traversal.pyx",
                "yt/utilities/lib/FixedInterpolator.c",
                "yt/utilities/lib/kdtree.c"],
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               extra_compile_args=omp_args,
               extra_link_args=omp_args,
               depends = ["yt/utilities/lib/VolumeIntegrator.pyx",
                          "yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/kdtree.h",
                          "yt/utilities/lib/FixedInterpolator.h",
                          "yt/utilities/lib/fixed_interpolator.pxd",
                          "yt/utilities/lib/field_interpolation_tables.pxd",
                          ]
          )
    config.add_extension("write_array",
                         ["yt/utilities/lib/write_array.pyx"])
    config.add_extension("GridTree", 
    ["yt/utilities/lib/GridTree.pyx"],
        libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("amr_kdtools", 
                         ["yt/utilities/lib/amr_kdtools.pyx"],
                         libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_subpackage("tests")

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
