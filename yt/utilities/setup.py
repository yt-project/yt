#!/usr/bin/env python


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('utilities', parent_package, top_path)
    config.add_subpackage("amr_kdtree")
    config.add_subpackage("poster")
    config.add_subpackage("answer_testing")
    config.add_subpackage("spatial")
    config.add_subpackage("grid_data_format")
    config.add_subpackage("parallel_tools")
    config.add_subpackage("lib")
    config.add_extension("data_point_utilities",
                         "yt/utilities/data_point_utilities.c",
                         libraries=["m"])
    config.add_subpackage("tests")
    config.add_subpackage("pyparselibconfig")
    config.make_config_py()  # installs __config__.py
    # config.make_svn_version_py()
    return config
