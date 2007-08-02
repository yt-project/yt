#!/usr/bin/env python
import os, sys, os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lagos',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.add_extension("PointCombine", "src/lagos/PointCombine.c", libraries=["m"])
    sys.argv.extend(["config_fc","--f77flags","'-Dr16 -ffixed-line-length-none -fno-second-underscore -DPYFORT -DNOMETALS'"])
    config.add_extension("EnzoFortranRoutines", \
                        ["src/lagos/solve_rate_cool.pyf", "src/lagos/f_src/*.F"], \
                        include_dirs=["/u/ki/mturk/Research/enzo/src/"])
    return config
