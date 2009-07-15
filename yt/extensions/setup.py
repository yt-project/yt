#!/usr/bin/env python
import setuptools

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('extensions',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.make_svn_version_py()
    config.add_subpackage("lightcone")
    config.add_subpackage("spatial")
    return config
