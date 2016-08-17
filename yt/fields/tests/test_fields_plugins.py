# -*- coding: UTF-8 -*-
#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import unittest
import yt
from yt.config import ytcfg, CONFIG_DIR
from yt.testing import fake_random_ds

TEST_PLUGIN_FILE = '''def _myfunc(field, data):
    return np.random.random(data['density'].shape)
add_field('random', dimensions='dimensionless',
          function=_myfunc, units='auto')'''


def setUpModule():
    my_plugin_name = ytcfg.get('yt', 'pluginfilename')
    # In the following order if pluginfilename is: an absolute path, located in
    # the CONFIG_DIR, located in an obsolete config dir.
    old_config_dir = os.path.join(os.path.expanduser('~'), '.yt')
    for base_prefix in ('', CONFIG_DIR, old_config_dir):
        potential_plugin_file = os.path.join(base_prefix, my_plugin_name)
        if os.path.isfile(potential_plugin_file):
            os.rename(potential_plugin_file,
                      potential_plugin_file + '.bak_test')

    plugin_file = os.path.join(CONFIG_DIR, my_plugin_name)
    with open(plugin_file, 'w') as fh:
        fh.write(TEST_PLUGIN_FILE)


def tearDownModule():
    my_plugin_name = ytcfg.get('yt', 'pluginfilename')
    plugin_file = os.path.join(CONFIG_DIR, my_plugin_name)
    os.remove(plugin_file)

    old_config_dir = os.path.join(os.path.expanduser('~'), '.yt')
    for base_prefix in ('', CONFIG_DIR, old_config_dir):
        potential_plugin_file = os.path.join(base_prefix, my_plugin_name)
        if os.path.isfile(potential_plugin_file + '.bak_test'):
            os.rename(potential_plugin_file + '.bak_test',
                      potential_plugin_file)


class TestPluginFile(unittest.TestCase):

    def testCustomField(self):
        plugin_file = os.path.join(
            CONFIG_DIR, ytcfg.get('yt', 'pluginfilename'))
        msg = 'INFO:yt:Loading plugins from %s' % plugin_file

        with self.assertLogs('yt', level='INFO') as cm:
            yt.enable_plugins()
            self.assertEqual(cm.output, [msg])

        ds = fake_random_ds(16)
        dd = ds.all_data()
        self.assertEqual(str(dd['random'].units), 'dimensionless')
        self.assertEqual(dd['random'].shape, dd['density'].shape)
