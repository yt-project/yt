#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('ramses', parent_package, top_path)
    config.add_extension("_ramses_reader",
        ["yt/frontends/ramses/_ramses_reader.pyx"],
        language="c++",
        include_dirs=["yt/frontends/ramses/ramses_headers/"],
        libraries=["stdc++"],
        depends=glob.glob("yt/frontends/ramses/ramses_headers/*.hh")
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
