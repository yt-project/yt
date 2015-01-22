#!/usr/bin/env python
import setuptools
import os, sys, os.path, glob, \
    tempfile, subprocess, shutil

def check_for_openmp():
    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    if os.name == 'nt': return False

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')

        # Attempt to compile a test script.
        # See http://openmp.org/wp/openmp-compilers/
        filename = r'test.c'
        file = open(filename,'wt', 1)
        file.write(
            "#include <omp.h>\n"
            "#include <stdio.h>\n"
            "int main() {\n"
            "#pragma omp parallel\n"
            "printf(\"Hello from thread %d, nthreads %d\\n\", omp_get_thread_num(), omp_get_num_threads());\n"
            "}"
            )
        file.flush()
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
                ["yt/utilities/lib/ContourFinding.pyx"],
                include_dirs=["yt/utilities/lib/",
                              "yt/geometry/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/utilities/lib/amr_kdtools.pxd",
                         "yt/utilities/lib/ContourFinding.pxd",
                         "yt/geometry/oct_container.pxd"])
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
    config.add_extension("origami", 
                ["yt/utilities/lib/origami.pyx",
                 "yt/utilities/lib/origami_tags.c"],
                include_dirs=["yt/utilities/lib/"],
                depends=["yt/utilities/lib/origami_tags.h"])
    config.add_extension("image_utilities", 
                         ["yt/utilities/lib/image_utilities.pyx"],
                         libraries=["m"],
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
               depends = ["yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/kdtree.h",
                          "yt/utilities/lib/FixedInterpolator.h",
                          "yt/utilities/lib/fixed_interpolator.pxd",
                          "yt/utilities/lib/field_interpolation_tables.pxd",
                          ]
          )
    config.add_extension("write_array",
                         ["yt/utilities/lib/write_array.pyx"])
    config.add_extension("ragged_arrays",
                         ["yt/utilities/lib/ragged_arrays.pyx"])
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
        print("INCLUDE AND LIB DIRS", idir, ldir)
        config.add_extension("perftools_wrap",
                ["yt/utilities/lib/perftools_wrap.pyx"],
                libraries=["profiler"],
                library_dirs = [ldir],
                include_dirs = [idir],
            )
    config.make_config_py() # installs __config__.py
    return config
