"""
This module is very simple.  It imports the configuration
we have written for yt.
Everything will be returned in a global config dictionary ``ytcfg``

"""

from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
from yt.extern.six.moves import configparser

ytcfg_defaults = dict(
    serialize = 'False',
    onlydeserialize = 'False',
    timefunctions = 'False',
    logfile = 'False',
    coloredlogs = 'False',
    suppressstreamlogging = 'False',
    stdoutStreamLogging = 'False',
    loglevel = '20',
    inline = 'False',
    numthreads = '-1',
    __withintesting = 'False',
    __parallel = 'False',
    __global_parallel_rank = '0',
    __global_parallel_size = '1',
    __topcomm_parallel_rank = '0',
    __topcomm_parallel_size = '1',
    __command_line = 'False',
    storeparameterfiles = 'False',
    parameterfilestore = 'parameter_files.csv',
    maximumstoreddatasets = '500',
    skip_dataset_cache = 'True',
    loadfieldplugins = 'True',
    pluginfilename = 'my_plugins.py',
    parallel_traceback = 'False',
    pasteboard_repo = '',
    reconstruct_index = 'False',
    test_storage_dir = '/does/not/exist',
    test_data_dir = '/does/not/exist',
    enzo_db = '',
    hub_url = 'https://hub.yt-project.org/upload',
    hub_api_key = '',
    notebook_password = '',
    answer_testing_tolerance = '3',
    answer_testing_bitwise = 'False',
    gold_standard_filename = 'gold311',
    local_standard_filename = 'local001',
    answer_tests_url = 'http://answers.yt-project.org/{1}_{2}',
    sketchfab_api_key = 'None',
    imagebin_api_key = 'e1977d9195fe39e',
    imagebin_upload_url = 'https://api.imgur.com/3/upload',
    imagebin_delete_url = 'https://api.imgur.com/3/image/{delete_hash}',
    thread_field_detection = 'False',
    ignore_invalid_unit_operation_errors = 'False',
    chunk_size = '1000',
    xray_data_dir = '/does/not/exist',
    default_colormap = 'arbre',
    ray_tracing_engine = 'embree',
    )
# Here is the upgrade.  We're actually going to parse the file in its entirety
# here.  Then, if it has any of the Forbidden Sections, it will be rewritten
# without them.

__fn = os.path.expanduser("~/.yt/config")
if os.path.exists(__fn):
    f = open(__fn).read()
    if any(header in f for header in ["[lagos]","[raven]","[fido]","[enki]"]):
        print("***********************************************************")
        print("* Upgrading configuration file to new format; saving old. *")
        print("***********************************************************")
        # This is of the old format
        cp = configparser.ConfigParser()
        cp.read(__fn)
        # NOTE: To avoid having the 'DEFAULT' section here,
        # we are not passing in ytcfg_defaults to the constructor.
        new_cp = configparser.ConfigParser()
        new_cp.add_section("yt")
        for section in cp.sections():
            for option in cp.options(section):
                # We changed them all to lowercase
                if option.lower() in ytcfg_defaults:
                    new_cp.set("yt", option, cp.get(section, option))
                    print("Setting %s to %s" % (option, cp.get(section, option)))
        open(__fn + ".old", "w").write(f)
        new_cp.write(open(__fn, "w"))
# Pathological check for Kraken
#elif os.path.exists("~/"):
#    if not os.path.exists("~/.yt"):
#            print "yt is creating a new directory, ~/.yt ."
#            os.mkdir(os.path.exists("~/.yt/"))
#    # Now we can read in and write out ...
#    new_cp = configparser.ConfigParser(ytcfg_defaults)
#    new_cp.write(__fn)

class YTConfigParser(configparser.ConfigParser):
    def __setitem__(self, key, val):
        self.set(key[0], key[1], val)
    def __getitem__(self, key):
        self.get(key[0], key[1])

if os.path.exists(os.path.expanduser("~/.yt/config")):
    ytcfg = YTConfigParser(ytcfg_defaults)
    ytcfg.read(['yt.cfg', os.path.expanduser('~/.yt/config')])
else:
    ytcfg = YTConfigParser(ytcfg_defaults)
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
