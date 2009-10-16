"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: ytcfg

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import ConfigParser, os, os.path, types


ytcfgDefaults = {
    "fido":{
        'RunDir': os.path.join(os.getenv("HOME"),'.yt/EnzoRuns/'),
        'WaitBetween':'5',
        'OtherFiles':'rates.out,cool_rates.out',
        'NewOutputCreated':'newOutput',
        'GlobPatterns':'*.hierarchy,*.dir',
        'NewDirectoryPattern':'%s.dir'
        },
    "reason":{
        'width':"600",
        'height':"600",
        'centeronmax':'False',
        'minpbar':'300',
        },
    "lagos":{
        'ReconstructHierarchy': 'True',
        'serialize' : 'True',
        'onlydeserialize' : 'False',
        'loadfieldplugins':'False',
        'pluginfilename':'yt_plugins.py',
        },
    "yt":{
        'time_functions': 'False',
        'LogFile': 'False',
        'LogFileName': 'yt.log',
        'suppressStreamLogging': 'False',
        'LogLevel': '20',
        'timefunctions':'False',
        'inGui':'False',
        'inline':'False',
        '__parallel':'False',
        '__parallel_rank':'0',
        '__parallel_size':'1',
        'StoreParameterFiles': 'True',
        'ParameterFileStore': 'parameter_files.csv',
        'MaximumStoredPFs': '500',
         },
    "raven":{
        'ImagePath':".",
        'ImageSkel': '%(bn)s_%(width)010i_%(unit)s',
        'colormap': 'jet',
        }
    }

class YTConfigParser(ConfigParser.ConfigParser):
    """
    Simple class providing some functionality I wish existed in the ConfigParser
    module already
    """
    def __init__(self, fn, defaults=None):
        if not defaults: defaults = {}
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
        ConfigParser.ConfigParser.set(self, section, opt, val)
    def __getitem__(self, item):
        if hasattr(item,'__getitem__'):
            tr = []
            for it in item[1:]:
                tr.append(self.get(item[0], it, raw=True))
            if len(tr) == 1:
                return tr[0]
            return tr
        else:
            raise KeyError
    def __setitem__(self, item, val):
        if not isinstance(item, types.TupleType) or not len(item) == 2:
            raise KeyError
        self.set(item[0], item[1], val)

if os.path.exists(os.path.expanduser("~/.yt/config")):
    ytcfg = YTConfigParser(['yt.cfg', os.path.expanduser('~/.yt/config')],
                           ytcfgDefaults)
else:
    ytcfg = YTConfigParser(['yt.cfg'],
                        ytcfgDefaults)

# Now we have parsed the config file.  Overrides come from the command line.

# This should be implemented at some point.  The idea would be to have a set of
# command line options, fed through an option parser, that would override
# the settings in ytcfg.  *However*, because we want to have the command-line
# scripts work, we'd probably want to have them only be long options, and also
# along the lines of --yt-something=somethingelse.  The command line scripts
# would then not get their options from sys.argv, but instead from this module.
