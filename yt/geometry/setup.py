#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('geometry',parent_package,top_path)
    config.add_extension("oct_container", 
                ["yt/geometry/oct_container.pyx"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.add_extension("selection_routines", 
                ["yt/geometry/selection_routines.pyx"],
                extra_compile_args=['-fopenmp'],
                extra_link_args=['-fopenmp'],
                include_dirs=["yt/utilities/lib/"],
                libraries=["m"], depends=["yt/utilities/lib/fp_utils.pxd"])
    config.make_config_py() # installs __config__.py
    #config.make_svn_version_py()
    return config
