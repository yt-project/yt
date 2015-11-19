#!/usr/bin/env python

def configuration(parent_package = '', top_path = None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('spatial', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_sconscript('SConstruct')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer = "SciPy Developers",
          author = "Anne Archibald",
          maintainer_email = "scipy-dev@scipy.org",
          description = "Spatial algorithms and data structures",
          url = "http://www.scipy.org",
          license = "SciPy License (BSD Style)",
          **configuration(top_path='').todict()
          )
