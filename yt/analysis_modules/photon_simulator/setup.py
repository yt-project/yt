#!/usr/bin/env python


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('photon_simulator', parent_package, top_path)
    config.add_extension("utils",
                         ["yt/analysis_modules/photon_simulator/utils.pyx"])
    config.add_subpackage("tests")
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
