#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('volume_rendering',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.make_svn_version_py()
    config.add_extension("VolumeIntegrator", "yt/extensions/volume_rendering/VolumeIntegrator.c")
    return config
