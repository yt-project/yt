#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('raven',parent_package,top_path)
    config.add_subpackage("deliveration")
    config.make_config_py() # installs __config__.py
    config.add_extension("_MPL", "_MPL.c", libraries=["m"])
    return config
