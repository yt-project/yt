#!/usr/bin/env python
import setuptools
import os, sys, os.path
import glob


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('artio', parent_package, top_path)

    config.add_extension("_artio_reader",
        ["yt/frontends/artio/_artio_reader.pyx"],
        language="c",
        include_dirs=["yt/frontends/artio/artio_headers/"],
        libraries=["stdc"], #snl? maybe m or something else? ramses = stdc++
        depends=glob.glob("yt/frontends/artio/artio_headers/*.hh") #snl don't know glob.glob or *.hh's
        )

    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
