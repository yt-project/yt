#!/usr/bin/env python
import os.path


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('halo_finding', parent_package, top_path)
    config.add_subpackage("fof")
    config.add_subpackage("hop")
    if os.path.exists("rockstar.cfg"):
        config.add_subpackage("rockstar")
    config.make_config_py()  # installs __config__.py
    return config
