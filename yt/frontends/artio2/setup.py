#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
from glob import glob

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('artio2', parent_package, top_path)
    config.add_extension("_artio_reader",
        ["yt/frontends/artio2/_artio_reader.pyx"]+
            glob("yt/frontends/artio2/artio/*.c"),
        language="c",
        include_dirs=["yt/frontends/artio2/artio/"],
        depends=glob("yt/frontends/artio2/artio/*.h") )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
