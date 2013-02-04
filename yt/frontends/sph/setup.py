#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('sph', parent_package, top_path)
    config.add_extension("smoothing_kernel",
        ["yt/frontends/sph/smoothing_kernel.pyx"],
        include_dirs=["yt/frontends/sph/",
                      "yt/geometry/",
                      "yt/utilities/lib/"],
        depends=glob.glob("yt/geometry/*.pxd"),
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
