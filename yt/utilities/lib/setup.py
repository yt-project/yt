#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import tempfile
import subprocess
import shutil
import pkg_resources
from sys import platform as _platform


def check_for_openmp():
    """Returns True if local setup supports OpenMP, False otherwise"""

    # See https://bugs.python.org/issue25150
    if sys.version_info[:3] == (3, 5, 0):
        return False

    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    exit_code = 1

    if os.name == 'nt': return False

    try:
        os.chdir(tmpdir)

        # Get compiler invocation
        compiler = os.getenv('CC', 'cc')
        compiler = compiler.split(' ')

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
            exit_code = subprocess.call(compiler + ['-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)

        # Clean up
        file.close()
    except OSError:
        return False
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

    return exit_code == 0

def check_for_pyembree():
    try:
        fn = pkg_resources.resource_filename("pyembree", "rtcore.pxd")
    except ImportError:
        return None
    return os.path.dirname(fn)


def read_embree_location():
    '''

    Attempts to locate the embree installation. First, we check for an
    EMBREE_DIR environment variable. If one is not defined, we look for
    an embree.cfg file in the root yt source directory. Finally, if that 
    is not present, we default to /usr/local. If embree is installed in a
    non-standard location and none of the above are set, the compile will
    not succeed. This only gets called if check_for_pyembree() returns 
    something other than None.

    '''

    rd = os.environ.get('EMBREE_DIR')
    if rd is not None:
        return rd
    print("EMBREE_DIR not set. Attempting to read embree.cfg")
    try:
        rd = open("embree.cfg").read().strip()
        return rd
    except IOError:
        print("Reading Embree location from embree.cfg failed.")
        print("If compilation fails, please place the base directory")
        print("of your Embree install in embree.cfg and restart.")
        return '/usr/local'

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lib',parent_package,top_path)
    if check_for_openmp() is True:
        omp_args = ['-fopenmp']
    else:
        omp_args = None
    if check_for_pyembree() is not None:
        loc = read_embree_location()
        embree_include_dir = loc + '/include'
        embree_lib_dir = loc + '/lib'
        embree_args = ['-I/' + embree_include_dir, '-L/' + embree_lib_dir]
        if _platform == "darwin":
            embree_lib_name = "embree.2"
        else:
            embree_lib_name = "embree"
    else:
        embree_args = None
    # Because setjmp.h is included by lots of things, and because libpng hasn't
    # always properly checked its header files (see
    # https://bugzilla.redhat.com/show_bug.cgi?id=494579 ) we simply disable
    # support for setjmp.
    config.add_extension("bitarray", 
                ["yt/utilities/lib/bitarray.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/bitarray.pxd"])
    config.add_extension("bounding_volume_hierarchy", 
                         ["yt/utilities/lib/bounding_volume_hierarchy.pyx"],
                         libraries=["m"], 
                         depends=["yt/utilities/lib/bounding_volume_hierarchy.pxd"])
    config.add_extension("particle_mesh_operations", 
                ["yt/utilities/lib/particle_mesh_operations.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("contour_finding", 
                ["yt/utilities/lib/contour_finding.pyx"],
                include_dirs=["yt/utilities/lib/",
                              "yt/geometry/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/utilities/lib/amr_kdtools.pxd",
                         "yt/utilities/lib/grid_traversal.pxd",
                         "yt/utilities/lib/contour_finding.pxd",
                         "yt/geometry/oct_container.pxd"])
    config.add_extension("depth_first_octree", 
                ["yt/utilities/lib/depth_first_octree.pyx"],
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
    config.add_extension("interpolators", 
                ["yt/utilities/lib/interpolators.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("alt_ray_tracers", 
                ["yt/utilities/lib/alt_ray_tracers.pyx"],
                libraries=["m"], depends=[])
    config.add_extension("marching_cubes", 
                ["yt/utilities/lib/marching_cubes.pyx",
                 "yt/utilities/lib/fixed_interpolator.c"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/utilities/lib/fixed_interpolator.pxd",
                         "yt/utilities/lib/fixed_interpolator.h",
                ])
    config.add_extension("misc_utilities", 
                ["yt/utilities/lib/misc_utilities.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("pixelization_routines",
                ["yt/utilities/lib/pixelization_routines.pyx",
                 "yt/utilities/lib/pixelization_constants.c"],
               include_dirs=["yt/utilities/lib/"],
               language="c++",
               libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd",
                                         "yt/utilities/lib/pixelization_constants.h",
                                         "yt/utilities/lib/element_mappings.pxd"])
    config.add_extension("basic_octree", 
                ["yt/utilities/lib/basic_octree.pyx"],
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
    config.add_extension("points_in_volume", 
                ["yt/utilities/lib/points_in_volume.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("quad_tree", 
                ["yt/utilities/lib/quad_tree.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("ray_integrators", 
                ["yt/utilities/lib/ray_integrators.pyx"],
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
                "yt/utilities/lib/fixed_interpolator.c",
                "yt/utilities/lib/kdtree.c"],
               include_dirs=["yt/utilities/lib/"],
               libraries=["m"], 
               extra_compile_args=omp_args,
               extra_link_args=omp_args,
               depends = ["yt/utilities/lib/fp_utils.pxd",
                          "yt/utilities/lib/kdtree.h",
                          "yt/utilities/lib/fixed_interpolator.h",
                          "yt/utilities/lib/fixed_interpolator.pxd",
                          "yt/utilities/lib/field_interpolation_tables.pxd",
                          ]
          )
    config.add_extension("write_array",
                         ["yt/utilities/lib/write_array.pyx"])
    config.add_extension("element_mappings",
                         ["yt/utilities/lib/element_mappings.pyx"],
                         libraries=["m"], depends=["yt/utilities/lib/element_mappings.pxd"])
    config.add_extension("ragged_arrays",
                         ["yt/utilities/lib/ragged_arrays.pyx"])
    config.add_extension("amr_kdtools", 
                         ["yt/utilities/lib/amr_kdtools.pyx"],
                         libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("line_integral_convolution",
                         ["yt/utilities/lib/line_integral_convolution.pyx"])

    include_dirs = check_for_pyembree()
    if include_dirs is not None:
        config.add_extension("mesh_construction",
                             ["yt/utilities/lib/mesh_construction.pyx"],
                             include_dirs=["yt/utilities/lib", include_dirs],
                             libraries=["m", embree_lib_name], language="c++",
                             extra_compile_args=embree_args,
                             extra_link_args=embree_args,
                             depends=["yt/utilities/lib/mesh_construction.pxd"])
        config.add_extension("mesh_traversal",
                             ["yt/utilities/lib/mesh_traversal.pyx"],
                             include_dirs=["yt/utilities/lib", include_dirs],
                             libraries=["m", embree_lib_name], language="c++",
                             extra_compile_args=embree_args,
                             extra_link_args=embree_args,
                             depends=["yt/utilities/lib/mesh_traversal.pxd",
                                      "yt/utilities/lib/grid_traversal.pxd"])
        config.add_extension("mesh_samplers",
                             ["yt/utilities/lib/mesh_samplers.pyx"],
                             include_dirs=["yt/utilities/lib", include_dirs],
                             libraries=["m", embree_lib_name], language="c++",
                             extra_compile_args=embree_args,
                             extra_link_args=embree_args,
                             depends=["yt/utilities/lib/mesh_samplers.pxd",
                                      "yt/utilities/lib/element_mappings.pxd",
                                      "yt/utilities/lib/mesh_construction.pxd"])
        config.add_extension("mesh_intersection",
                             ["yt/utilities/lib/mesh_intersection.pyx"],
                             include_dirs=["yt/utilities/lib", include_dirs],
                             libraries=["m", embree_lib_name], language="c++",
                             extra_compile_args=embree_args,
                             extra_link_args=embree_args,
                             depends=["yt/utilities/lib/mesh_intersection.pxd",
                                      "yt/utilities/lib/mesh_construction.pxd"])
    config.add_subpackage("tests")

    if os.environ.get("GPERFTOOLS", "no").upper() != "NO":
        gpd = os.environ["GPERFTOOLS"]
        idir = os.path.join(gpd, "include")
        ldir = os.path.join(gpd, "lib")
        print(("INCLUDE AND LIB DIRS", idir, ldir))
        config.add_extension("perftools_wrap",
                ["yt/utilities/lib/perftools_wrap.pyx"],
                libraries=["profiler"],
                library_dirs = [ldir],
                include_dirs = [idir],
            )
    config.make_config_py() # installs __config__.py
    return config
