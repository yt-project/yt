import os, os.path
import sys
import time
import subprocess
import ez_setup
ez_setup.use_setuptools()

import setuptools

APP = ['reason.py']
DATA_FILES = []
PY2APP_OPTIONS = {'argv_emulation': True}
VERSION = "1.5dev"

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
        version = VERSION,
        description = "An analysis and visualization toolkit for Adaptive Mesh " \
                    + "Refinement data, specifically for the Enzo and Orion codes.",
        classifiers = [],
        keywords='',
        install_requires = ['matplotlib', 'numpy','ipython'],
        extras_require = { 'GUI' : ['wxPython'],
                           'storage' : ['tables'], },
        author="Matthew J. Turk",
        author_email="matthewturk@gmail.com",
        url = "http://yt.enzotools.org/",
        license="GPL-3",
        configuration=configuration,
        app=APP, # for py2app
        data_files=DATA_FILES,
        options={'py2app':PY2APP_OPTIONS},
        zip_safe=False,
        )
    return

if __name__ == '__main__':
    setup_package()
