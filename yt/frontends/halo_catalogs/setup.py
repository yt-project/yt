#!/usr/bin/env python
from numpy.distutils.misc_util import Configuration


def configuration(parent_package='', top_path=None):
    config = Configuration('halo_catalogs', parent_package, top_path)
    config.add_subpackage("halo_catalog")
    config.add_subpackage("owls_subfind")
    config.add_subpackage("rockstar")
    config.make_config_py()
    return config
