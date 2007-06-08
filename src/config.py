"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: ytcfg

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Implement default rc file outputting
"""

import ConfigParser, os, os.path, sys

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
        'LogLevel': '30', \
        'unifiedlogfile': '1', \
        'User':os.getenv("USER"), \
         }, \
    "raven":{
        'ImagePath':".", \
        'ImageSkel': '%(bn)s_%(width)010i_%(unit)s'\
        } \
    }

class YTConfigParser(ConfigParser.ConfigParser):
    def __init__(self, fn, defaults={}):
        ConfigParser.ConfigParser.__init__(self)
        # Note that we're not going to pass in defaults
        self.read(fn)
        # Okay, we're populated.  Now, we will in the additional values
        # as needed.
        for section in defaults.keys():
            opts = defaults[section]
            if not self.has_section(section):
                self.add_section(section)
            for opt, val in opts.items():
                if not self.has_option(section, opt):
                    self.set(section, opt, val)

#ytcfg = ConfigParser.ConfigParser(defaults)
ytcfg = YTConfigParser(['yt.cfg', os.path.expanduser('~/.yt/config')], ytcfgDefaults)
#ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])

# I am hereby getting rid of usermodules...  (They used to be below.)  I can
# actually think of no reasonable purpose that they would serve.  Also, it
# messed up the logging import.
