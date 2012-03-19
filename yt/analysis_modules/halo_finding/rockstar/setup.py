#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('rockstar',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    #config.make_svn_version_py()
    rd = os.environ["ROCKSTAR_DIR"]
    config.add_extension("rockstar_interface",
                         "yt/analysis_modules/halo_finding/rockstar/rockstar_interface.pyx",
                         library_dirs=[rd],
                         libraries=["rockstar"],
                         include_dirs=[rd,
                                       os.path.join(rd, "io"),
                                       os.path.join(rd, "util")])
    return config

