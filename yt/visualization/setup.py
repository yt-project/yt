#!/usr/bin/env python
import setuptools


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('visualization', parent_package, top_path)
    config.add_subpackage("image_panner")
    config.add_subpackage("volume_rendering")
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    config.add_extension("_MPL", "_MPL.c", libraries=["m"])
    return config
