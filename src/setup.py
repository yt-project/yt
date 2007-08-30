#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('yt',parent_package,top_path)
    config.add_subpackage('lagos')
    config.add_subpackage('raven')
    config.add_subpackage('enki')
    config.add_subpackage('fido')
    config.add_subpackage('reason')
    config.make_config_py()
    return config

if __name__ == '__main__':
    # Remove current working directory from sys.path
    # to avoid importing numpy.distutils as Python std. distutils:
    import os, sys
    from numpy.distutils.core import setup
    setup(configuration=configuration)
