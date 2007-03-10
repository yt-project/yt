"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: YTConfig

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import ConfigParser, os, os.path

defaults = {'RunDir': os.path.join(os.getenv("HOME"),'EnzoRuns'),\
            'WorkDir': os.path.join('/usr/work/', os.getenv("USER")),\
            'EnzoInterfacePath': '/usr/work/mturk/local/lib/python2.5/site-packages'}

ytcfg = ConfigParser.ConfigParser(defaults)
ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])
