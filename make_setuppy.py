ss = r"""#!/usr/bin/env python
import setuptools
import os, sys, os.path

import os.path

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('%s',parent_package,top_path)
    config.make_config_py() # installs __config__.py
    config.make_svn_version_py()
    return config
"""

import sys, os
n = sys.argv[-1]

fn = os.path.join(os.path.abspath(n), "setup.py")
if os.path.exists(fn):
    print "%s exists!  Not overwriting." % fn
bn = os.path.basename(n)

print "%s -> %s" % (bn, fn)

open(fn, "w").write(ss % bn)
