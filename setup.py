#from numpy.distutils.core import setup, Extension, Configuration

#_ec = Extension('yt.lagos.EnzoCombine',['src/lagos/EnzoCombine.c'], \
                #include_dirs=[numpyincludedirs], \
                #libraries=['m'])
#DOCLINES = __doc__.split("\n")

import os
import sys
import time

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    
    config.make_config_py()
    config.add_subpackage('yt','src')

    #config.add_data_files(('numpy',['*.txt','COMPATIBILITY',
                                    #'scipy_compatibility']))

    #config.get_version('numpy/version.py') # sets config.version
    
    return config

def setup_package():

    from numpy.distutils.core import setup

    setup(
        name = "yt",
        version = time.strftime("%y%m%d"),
        description = "A set of classes for manipulating Enzo data",
        url = "http://www.stanford.edu/~mturk/raven.html",
        author="Matthew Turk",
        author_email="mturk@stanford.edu",
        license="GPL-2",
        configuration=configuration
        )
    return

if __name__ == '__main__':
    setup_package()

        #package_dir={"yt":"src"},
        #packages = ["yt.enki", "yt.enki.mes", "yt.fido", "yt.lagos", "yt", "yt.raven", "yt.deliverator"],
        #scripts = ["scripts/fdigup","scripts/fimport","scripts/frevert","scripts/fbranch", 
                   #"scripts/ffetch", "scripts/yt", "scripts/fbury","scripts/fido"],
        #data_files = [('share/doc/yt/', ['examples/test_enki.py',
          #'examples/test_fido.py', 'examples/test_raven.py', 'examples/test_lagos.py',
          #'examples/get_particles_enzorun.py'])],
