"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: YTConfig

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import ConfigParser, os

ytcfg = ConfigParser.ConfigParser()
ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])
