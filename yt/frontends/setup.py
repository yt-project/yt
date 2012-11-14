#!/usr/bin/env python
import setuptools


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('frontends', parent_package, top_path)
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    config.add_subpackage("athena")
    config.add_subpackage("gdf")
    config.add_subpackage("chombo")
    config.add_subpackage("enzo")
    config.add_subpackage("flash")
    config.add_subpackage("nyx")
    config.add_subpackage("orion")
    config.add_subpackage("artio")
    config.add_subpackage("artio2")
    config.add_subpackage("ramses")
    config.add_subpackage("tiger")
    config.add_subpackage("art")
    config.add_subpackage("maestro")
    config.add_subpackage("castro")
    config.add_subpackage("stream")
    return config
