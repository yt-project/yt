#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('geometry',parent_package,top_path)
    config.add_subpackage('coordinates')
    config.add_extension("grid_visitors", 
                ["yt/geometry/grid_visitors.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/grid_visitors.pxd"])
    config.add_extension("grid_container", 
                ["yt/geometry/grid_container.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/grid_container.pxd",
                         "yt/geometry/grid_visitors.pxd"])
    config.add_extension("oct_container", 
                ["yt/geometry/oct_container.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd"])
    config.add_extension("oct_visitors", 
                ["yt/geometry/oct_visitors.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd"])
    config.add_extension("particle_oct_container", 
                ["yt/geometry/particle_oct_container.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd"])
    config.add_extension("selection_routines", 
                ["yt/geometry/selection_routines.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/oct_visitors.pxd",
                         "yt/geometry/grid_container.pxd",
                         "yt/geometry/grid_visitors.pxd",
                         "yt/geometry/selection_routines.pxd"])
    config.add_extension("particle_deposit", 
                ["yt/geometry/particle_deposit.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd",
                         "yt/geometry/particle_deposit.pxd"])
    config.add_extension("particle_smooth", 
                ["yt/geometry/particle_smooth.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd",
                         "yt/geometry/particle_deposit.pxd",
                         "yt/geometry/particle_smooth.pxd"])
    config.add_extension("fake_octree", 
                ["yt/geometry/fake_octree.pyx"],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"],
                depends=["yt/utilities/lib/fp_utils.pxd",
                         "yt/geometry/oct_container.pxd",
                         "yt/geometry/selection_routines.pxd"])
    config.make_config_py() # installs __config__.py
    #config.make_svn_version_py()
    return config
