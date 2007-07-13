"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: ytcfg

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Implement default rc file outputting
"""

import ConfigParser, os, os.path, sys, types
from funcs import *

ytcfgDefaults = {
    "Fido":{
        'RunDir': os.path.join(os.getenv("HOME"),'EnzoRuns'),\
        'WorkDir': os.path.join('/usr/work/', os.getenv("USER")) \
        },\
    "SWIG":{
        'EnzoInterfacePath': '/usr/work/mturk/local/lib/python2.5/site-packages', \
        },\
    "yt":{
        'LogFile': '1', \
        'LogLevel': '20', \
        'unifiedlogfile': '1', \
        'User':os.getenv("USER"), \
        'timefunctions':'False'
         }, \
    "raven":{
        'ImagePath':".", \
        'ImageSkel': '%(bn)s_%(width)010i_%(unit)s',\
        'backend': 'HD'\
        } \
    }

class YTConfigParser(ConfigParser.ConfigParser):
    def __init__(self, fn, defaults={}):
        ConfigParser.ConfigParser.__init__(self)
        # Note that we're not going to pass in defaults
        self.read(fn)
        # Okay, we're populated.  Now, we will insert additional values
        # as needed.
        for section, opts in defaults.items():
            for opt, val in opts.items():
                if not self.has_option(section, opt):
                    self.set(section, opt, val)
    def set(self, section, opt, val):
        if not self.has_section(section):
            self.add_section(section)
        ConfigParser.ConfigParser.set(self,section, opt,val)
    def __getitem__(self, item):
        #if iterable(item):
        if hasattr(item,'__getitem__'):
            tr = []
            for it in item[1:]:
                tr.append(self.get(item[0],it,raw=True))
            if len(tr) == 1:
                return tr[0]
            return tr
        else:
            raise KeyError
    def __setitem__(self, item, val):
        if not isinstance(item, types.TupleType) or not len(item) == 2:
            raise KeyError
        self.set(item[0],item[1],val)

ytcfg = YTConfigParser(['yt.cfg', os.path.expanduser('~/.yt/config')], ytcfgDefaults)
