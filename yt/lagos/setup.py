#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def check_for_hdf5():
    if "HDF5_DIR" in os.environ:
        return os.environ["HDF5_DIR"]
    elif os.path.exists("hdf5.cfg"):
        return open("hdf5.cfg").read().strip().rstrip()
    print "Reading HDF5 location from hdf5.cfg failed."
    print "Please place the base directory of your HDF5 install in hdf5.cfg and restart."
    print "(ex: \"echo '/usr/local/' > hdf5.cfg\" )"
    sys.exit(1)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lagos',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.add_extension("PointCombine", "yt/lagos/PointCombine.c", libraries=["m"])
    #config.add_extension("RTIntegrator", "yt/lagos/RTIntegrator.c")
    config.add_extension("Interpolators", "yt/lagos/Interpolators.c")
    #config.add_extension("DepthFirstOctree", "yt/lagos/DepthFirstOctree.c")
    config.add_subpackage("hop")
    H5dir = check_for_hdf5()
    if H5dir is not None:
        include_dirs=[os.path.join(H5dir,"include")]
        library_dirs=[os.path.join(H5dir,"lib")]
        config.add_extension("HDF5LightReader", "yt/lagos/HDF5LightReader.c",
                             libraries=["m","hdf5"],
                             library_dirs=library_dirs, include_dirs=include_dirs)
    # Uncomment the next two lines if you want particle_density support
    #config.add_extension("cic_deposit", ["yt/lagos/enzo_routines/cic_deposit.pyf",
    #                                     "yt/lagos/enzo_routines/cic_deposit.f"])
    if 0:
        sys.argv.extend(["config_fc","--f77flags",
                         "'-Dr16 -ffixed-line-length-132 -fno-second-underscore -DPYFORT -DNOMETALS -ggdb -O0'"])
        config.add_extension("EnzoFortranRoutines", \
                            ["yt/lagos/solve_rate_cool.pyf", "yt/lagos/f_src/*.F"])
    return config
