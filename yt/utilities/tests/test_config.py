# -*- coding: UTF-8 -*-
#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import contextlib
import os
import sys
if sys.version_info.major < 3:
    try:
        import mock
    except ImportError:
        mock = None
else:
    import unittest.mock as mock
import unittest
import yt.utilities.command_line
import yt.config
from yt.config import \
    CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE, CONFIG_DIR, YTConfigParser
from yt.extern.six import StringIO
from yt.extern.six.moves.configparser import NoOptionError
from yt.fields.tests.test_fields_plugins import TEST_PLUGIN_FILE

_TEST_PLUGIN = '_test_plugin.py'
_DUMMY_CFG = ['[yt]', 'loglevel = 49', 'pluginfilename = ' + _TEST_PLUGIN]


@contextlib.contextmanager
def captureOutput():
    oldout, olderr = sys.stdout, sys.stderr
    try:
        out = [StringIO(), StringIO()]
        sys.stdout, sys.stderr = out
        yield out
    finally:
        sys.stdout, sys.stderr = oldout, olderr
        out[0] = out[0].getvalue()
        out[1] = out[1].getvalue()


class SysExitException(Exception):
    pass


def setUpModule():
    for cfgfile in (CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE):
        if os.path.exists(cfgfile):
            os.rename(cfgfile, cfgfile + '.bak_test')

            if cfgfile == CURRENT_CONFIG_FILE:
                yt.utilities.configure.CONFIG = YTConfigParser()
                if not yt.utilities.configure.CONFIG.has_section('yt'):
                    yt.utilities.configure.CONFIG.add_section('yt')


def tearDownModule():
    for cfgfile in (CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE): 
        if os.path.exists(cfgfile + '.bak_test'):
            os.rename(cfgfile + '.bak_test', cfgfile)


class TestYTConfig(unittest.TestCase):
    def _runYTConfig(self, args):
        args = ['yt', 'config'] + args
        retcode = 0

        with mock.patch.object(sys, 'argv', args),\
                mock.patch('sys.exit', side_effect=SysExitException) as exit,\
                captureOutput() as output:
            try:
                yt.utilities.command_line.run_main()
            except SysExitException:
                args = exit.mock_calls[0][1]
                retcode = args[0] if len(args) else 0
        return {
            'rc': retcode,
            'stdout': output[0],
            'stderr': output[1]
        }

    def _testKeyValue(self, key, val_set, val_get):
        info = self._runYTConfig(['set', 'yt', key, val_set])
        self.assertEqual(info['rc'], 0)

        info = self._runYTConfig(['get', 'yt', key])
        self.assertEqual(info['rc'], 0)
        self.assertEqual(info['stdout'].strip(), val_get)

        info = self._runYTConfig(['rm', 'yt', key])
        self.assertEqual(info['rc'], 0)

class TestYTConfigCommands(TestYTConfig):
    def testConfigCommands(self):
        # stub out test if mock isn't installed in Python2
        if mock is None:
            return

        self.assertFalse(os.path.exists(CURRENT_CONFIG_FILE))

        info = self._runYTConfig(['--help'])
        self.assertEqual(info['rc'], 0)
        self.assertEqual(info['stderr'], '')
        self.assertIn('Get and set configuration values for yt',
                      info['stdout'])

        info = self._runYTConfig(['list'])
        self.assertEqual(info['rc'], 0)
        self.assertIn('[yt]', info['stdout'])

        self._testKeyValue('__parallel', 'True', 'True')
        self._testKeyValue('test_data_dir', '~/yt-data',
                           os.path.expanduser('~/yt-data'))
        self._testKeyValue('test_data_dir', '$HOME/yt-data',
                           os.path.expandvars('$HOME/yt-data'))

        with self.assertRaises(NoOptionError):
            self._runYTConfig(['get', 'yt', 'foo'])
    
    def tearDown(self):
        if os.path.exists(CURRENT_CONFIG_FILE):
            os.remove(CURRENT_CONFIG_FILE)

class TestYTConfigMigration(TestYTConfig):

    def setUp(self):
        if not os.path.exists(os.path.dirname(_OLD_CONFIG_FILE)):
            os.makedirs(os.path.dirname(_OLD_CONFIG_FILE))
        
        with open(_OLD_CONFIG_FILE, 'w') as fh:
            for line in _DUMMY_CFG:
                fh.write('{}\n'.format(line))
        
        if os.path.exists(CURRENT_CONFIG_FILE):
            os.remove(CURRENT_CONFIG_FILE)
    
        my_plugin_name = _TEST_PLUGIN
        old_config_dir = os.path.join(os.path.expanduser('~'), '.yt')
        for base_prefix in ('', CONFIG_DIR, old_config_dir):
            potential_plugin_file = os.path.join(base_prefix, my_plugin_name)
            if os.path.isfile(potential_plugin_file):
                os.rename(potential_plugin_file,
                          potential_plugin_file + '.bak_test')

        plugin_file = os.path.join(old_config_dir, my_plugin_name)
        with open(plugin_file, 'w') as fh:
            fh.write(TEST_PLUGIN_FILE)

    def tearDown(self):
        if os.path.exists(CURRENT_CONFIG_FILE):
            os.remove(CURRENT_CONFIG_FILE)
        if os.path.exists(_OLD_CONFIG_FILE + '.bak'):
            os.remove(_OLD_CONFIG_FILE + '.bak')
        
        my_plugin_name = _TEST_PLUGIN
        plugin_file = os.path.join(CONFIG_DIR, my_plugin_name)
        if os.path.exists(plugin_file):
            os.remove(plugin_file)

        old_config_dir = os.path.join(os.path.expanduser('~'), '.yt')
        for base_prefix in ('', CONFIG_DIR, old_config_dir):
            potential_plugin_file = os.path.join(base_prefix, my_plugin_name)
            if os.path.isfile(potential_plugin_file + '.bak_test'):
                os.rename(potential_plugin_file + '.bak_test',
                          potential_plugin_file)

    def testConfigMigration(self):
        old_config_dir = os.path.dirname(_OLD_CONFIG_FILE)
        new_config_dir = os.path.dirname(CURRENT_CONFIG_FILE)

        self.assertFalse(os.path.exists(CURRENT_CONFIG_FILE))
        self.assertTrue(os.path.exists(_OLD_CONFIG_FILE))
        self.assertTrue(
            os.path.exists(os.path.join(old_config_dir, _TEST_PLUGIN)))
        
        info = self._runYTConfig(['migrate'])
        self.assertEqual(info['rc'], 0)

        self.assertTrue(os.path.exists(CURRENT_CONFIG_FILE))
        self.assertFalse(os.path.exists(_OLD_CONFIG_FILE))
        self.assertTrue(os.path.exists(_OLD_CONFIG_FILE + '.bak'))
        self.assertTrue(
            os.path.exists(os.path.join(new_config_dir, _TEST_PLUGIN)))
        self.assertTrue(
            os.path.exists(os.path.join(old_config_dir, _TEST_PLUGIN + '.bak'))
        )

        with open(CURRENT_CONFIG_FILE, 'r') as fh:
            new_cfg = ''.join(fh.readlines())
        self.assertEqual(new_cfg.strip().split('\n'), _DUMMY_CFG)
