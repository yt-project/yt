#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('hop',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.add_extension("EnzoHop", sources=
                                    ["EnzoHop.c",
                                     "hop_hop.c",
                                     "hop_kd.c",
                                     "hop_regroup.c",
                                     "hop_slice.c",
                                     "hop_smooth.c",])
    return config
