#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob

sourcefiles = ['artio_caller.pyx','artio_headers/artio_grid.c','artio_headers/artio_endian.c','artio_headers/artio_mpi.c', 'artio_headers/artio_parameter.c','artio_headers/artio_particle.c','artio_headers/artio_posix.c']

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('artio', parent_package, top_path)
    config.add_extension("_artio_reader",
        ["yt/frontends/artio/_artio_reader.pyx"],
        language="c++",
        include_dirs=["yt/frontends/artio/artio_headers/"],
        libraries=["stdc++"],
#        depends=glob.glob("yt/frontends/artio/artio_headers/*.hh")
        depends=glob.glob("yt/frontends/artio/artio_headers/*.hh,yt/frontends/artio/artio_headers/*.c")
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
