#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    
    artio_sources = glob.glob("yt/frontends/artio/artio_headers/*.c")

    config = Configuration('artio', parent_package, top_path)
    config.add_extension("_artio_caller",
        ["yt/frontends/artio/_artio_caller.pyx"]+artio_sources,
        include_dirs=["yt/frontends/artio/artio_headers/"],
        depends=artio_sources
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
