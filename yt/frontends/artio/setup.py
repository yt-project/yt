#!/usr/bin/env python
import setuptools
import os
import sys
import os.path
import glob

sourcefiles = ['yt/frontends/artio/_artio_caller.pyx',
               'yt/frontends/artio/artio_headers/artio_grid.c',
               'yt/frontends/artio/artio_headers/artio_endian.c',
               'yt/frontends/artio/artio_headers/artio_mpi.c', 
               'yt/frontends/artio/artio_headers/artio_parameter.c',
               'yt/frontends/artio/artio_headers/artio_particle.c',
               'yt/frontends/artio/artio_headers/artio_posix.c',
               'yt/frontends/artio/artio_headers/sfc.c',
               'yt/frontends/artio/artio_headers/sfc.h'
               ]

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('artio', parent_package, top_path)
    config.add_extension("_artio_caller",
        sourcefiles,
        include_dirs=["yt/frontends/artio/artio_headers/"],
        depends=glob.glob("yt/frontends/artio/artio_headers/*.c")
        )
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
