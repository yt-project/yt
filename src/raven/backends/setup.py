#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('backends',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.add_extension("_MPL", "src/raven/backends/_MPL.c", libraries=["m"])
    return config
