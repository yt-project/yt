#!/usr/bin/env python
import setuptools
import os
import sys


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('yt', parent_package, top_path)
    config.add_subpackage('analysis_modules')
    config.add_subpackage('astro_objects')
    config.add_subpackage('data_objects')
    config.add_subpackage('frontends')
    config.add_subpackage('gui')
    config.add_subpackage('utilities')
    config.add_subpackage('visualization')
    config.make_config_py()
    #config.make_svn_version_py()
    return config
