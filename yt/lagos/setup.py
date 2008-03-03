#!/usr/bin/env python
import os, sys, os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lagos',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.add_extension("PointCombine", "yt/lagos/PointCombine.c", libraries=["m"])
    if not os.environ.has_key("HDF5_INSTALL"):
        print "You need to set the environment variable HDF5_INSTALL to point to your HDF5 base"
        #raise KeyError("HDF5_INSTALL")
        sys.exit(1)
    H5dir = os.environ["HDF5_INSTALL"]
    include_dirs=[os.path.join(H5dir,"include")]
    library_dirs=[os.path.join(H5dir,"lib")]
    config.add_extension("HDF5LightReader", "yt/lagos/HDF5LightReader.c",
                         libraries=["m","hdf5_hl","hdf5"],
                         library_dirs=library_dirs, include_dirs=include_dirs)
    sys.argv.extend(["config_fc","--f77flags","'-Dr16 -ffixed-line-length-none -fno-second-underscore -DPYFORT -DNOMETALS -ggdb -O0'"])
    if 0:
        config.add_extension("EnzoFortranRoutines", \
                            ["yt/lagos/solve_rate_cool.pyf", "yt/lagos/f_src/*.F"])
    return config
