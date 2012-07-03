"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary: ytcfg

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

ytcfgDefaults = dict(
    serialize = 'True',
    onlydeserialize = 'False',
    timefunctions = 'False',
    logfile = 'False',
    coloredlogs = 'False',
    suppressstreamlogging = 'False',
    loglevel = '20',
    inline = 'False',
    numthreads = '-1',
    __withinreason = 'False',
    __parallel = 'False',
    __global_parallel_rank = '0',
    __global_parallel_size = '1',
    __topcomm_parallel_rank = '0',
    __topcomm_parallel_size = '1',
    __command_line = 'False',
    storeparameterfiles = 'False',
    parameterfilestore = 'parameter_files.csv',
    maximumstoredpfs = '500',
    loadfieldplugins = 'True',
    pluginfilename = 'my_plugins.py',
    parallel_traceback = 'False',
    pasteboard_repo = '',
    test_storage_dir = '/does/not/exist',
    enzo_db = '',
    hub_url = 'https://data.yt-project.org/upload',
    hub_api_key = '',
    )
# Here is the upgrade.  We're actually going to parse the file in its entirety
# here.  Then, if it has any of the Forbidden Sections, it will be rewritten
# without them.

__fn = os.path.expanduser("~/.yt/config")
if os.path.exists(__fn):
    f = open(__fn).read()
    if any(header in f for header in ["[lagos]","[raven]","[fido]","[enki]"]):
        print "***********************************************************"
        print "* Upgrading configuration file to new format; saving old. *"
        print "***********************************************************"
        # This is of the old format
        cp = ConfigParser.ConfigParser()
        cp.read(__fn)
        # NOTE: To avoid having the 'DEFAULT' section here,
        # we are not passing in ytcfgDefaults to the constructor.
        new_cp = ConfigParser.ConfigParser()
        new_cp.add_section("yt")
        for section in cp.sections():
            for option in cp.options(section):
                # We changed them all to lowercase
                if option.lower() in ytcfgDefaults:
                    new_cp.set("yt", option, cp.get(section, option))
                    print "Setting %s to %s" % (option, cp.get(section, option))
        open(__fn + ".old", "w").write(f)
        new_cp.write(open(__fn, "w"))
# Pathological check for Kraken
#elif os.path.exists("~/"):
#    if not os.path.exists("~/.yt"):
#            print "yt is creating a new directory, ~/.yt ."
#            os.mkdir(os.path.exists("~/.yt/"))
#    # Now we can read in and write out ...
#    new_cp = Configparser.ConfigParser(ytcfgDefaults)
#    new_cp.write(__fn)

class YTConfigParser(ConfigParser.ConfigParser):
    def __setitem__(self, key, val):
        self.set(key[0], key[1], val)

if os.path.exists(os.path.expanduser("~/.yt/config")):
    ytcfg = YTConfigParser(ytcfgDefaults)
    ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])
else:
    ytcfg = YTConfigParser(ytcfgDefaults)
    ytcfg.read(['yt.cfg'])
if not ytcfg.has_section("yt"):
    ytcfg.add_section("yt")

# Now we have parsed the config file.  Overrides come from the command line.

# This should be implemented at some point.  The idea would be to have a set of
# command line options, fed through an option parser, that would override
# the settings in ytcfg.  *However*, because we want to have the command-line
# scripts work, we'd probably want to have them only be long options, and also
# along the lines of --yt-something=somethingelse.  The command line scripts
# would then not get their options from sys.argv, but instead from this module.
