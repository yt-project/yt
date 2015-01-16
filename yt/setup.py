#!/usr/bin/env python
import setuptools
import os
import sys


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('yt', parent_package, top_path)
    if sys.version < '3':
        config.add_subpackage('analysis_modules') 
    else:
        print("no analysis modules with py3")

    config.add_subpackage('data_objects')
    config.add_subpackage('fields')
    config.add_subpackage('extern')
    config.add_subpackage('frontends')
    config.add_subpackage('geometry')
    config.add_subpackage('gui')
    config.add_subpackage('units')
    config.add_subpackage('utilities')
    config.add_subpackage('visualization')
    config.make_config_py()
    #config.make_svn_version_py()
    return config
