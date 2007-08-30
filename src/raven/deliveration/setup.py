#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('deliveration',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    return config
