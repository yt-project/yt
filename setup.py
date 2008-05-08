import os, os.path
import sys
import time
import subprocess
import ez_setup
ez_setup.use_setuptools()

import setuptools

APP = ['reason.py']
DATA_FILES = []
OPTIONS = {'argv_emulation': True}

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    
    config.make_config_py()
    config.add_subpackage('yt','yt')
    config.add_scripts("scripts/*")

    return config

def setup_package():

    from numpy.distutils.core import setup

    setup(
        name = "yt",
        version = "0.3",
        description = "A set of classes for manipulating Enzo Adaptive Mesh Refinement data",
        install_requires = ['matplotlib>=0.90.1',
                            'numpy>=1.0.3',
                            'wxPython>=2.8.7.1'],
        url = "http://yt.spacepope.org/",
        author="Matthew Turk",
        author_email="matt@yt.spacepope.org",
        license="GPL-3",
        configuration=configuration,
        #app=APP, data_files=DATA_FILES, options={'py2app':OPTIONS},
        )
    return

if __name__ == '__main__':
    setup_package()
