import os, os.path
import sys
import time
import subprocess

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
    config.add_scripts("scripts/*")

    return config

def setup_package():

    from numpy.distutils.core import setup

    setup(
        name = "yt",
        version = getVersion(),
        description = "A set of classes for manipulating Enzo data",
        url = "http://www.stanford.edu/~mturk/raven.html",
        author="Matthew Turk",
        author_email="mturk@stanford.edu",
        license="GPL-3",
        configuration=configuration
        )
    return

def getVersion():
    # First we try to SVN it
    repo = os.path.basename(os.path.realpath(__file__))
    p=subprocess.Popen("svn info %s" % repo, \
                              stdout=subprocess.PIPE, \
                              stderr=subprocess.STDOUT, \
                              shell=True, executable="/bin/bash")
    p.wait()
    version = None
    for line in p.stdout:
        if line.startswith("Revision:"):
            print "VERSION",line.split(":")[-1]
            version = "SVN-r%s" % (line.split(":")[-1].strip().rstrip())
    if not version: version = time.strftime("%y%m%d")
    return version

if __name__ == '__main__':
    setup_package()
