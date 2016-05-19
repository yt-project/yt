"""
A mechanism for listing available analysis modules.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import sys

def get_available_modules():
    modpath = os.path.abspath(os.path.dirname(__file__))
    available_modules = []
    for d in [os.path.join(modpath, f) for f in os.listdir(modpath)]:
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, "api.py")):
            available_modules.append(os.path.basename(d))
    return available_modules

class AnalysisModuleLoader(object):

    @property
    def available_modules(self):
        return get_available_modules()

    def __getattr__(self, attr):
        try:
            name = "yt.analysis_modules.%s.api" % (attr)
            __import__(name, level=-1)
            setattr(self, attr, sys.modules[name])
        except ImportError:
            raise AttributeError(attr)
        return getattr(self, attr)

    def __dir__(self):
        # This is a badly behaving object.  I was unable to get this line:
        #return super(AnalysisModuleLoader, self).__dir__() + self.available_modules
        # to work, so we simply return only the methods we know about.
        return ["available_modules"] + self.available_modules

amods = AnalysisModuleLoader()
