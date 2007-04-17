"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: YTConfig

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import ConfigParser, os, os.path, sys

defaults = {'RunDir': os.path.join(os.getenv("HOME"),'EnzoRuns'),\
            'WorkDir': os.path.join('/usr/work/', os.getenv("USER")),\
            'EnzoInterfacePath': '/usr/work/mturk/local/lib/python2.5/site-packages'}

ytcfg = ConfigParser.ConfigParser(defaults)
ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])

if ytcfg.has_section("usermodules") and ytcfg.has_option("usermodules","path"):
    umod_path = ytcfg.get("usermodules","path")
    mylog.info("Using %s as path to user modules", umod_path)
    sys.path = sys.path[:1] + [umod_path] + sys.path[1:] # We want '' to be the first
    # This is dangerous.  But, I know everyone is clamoring for it.
    for key in ytcfg.options("usermodules"):
        exec("import %s" % ytcfg.get("usermodules",key))
