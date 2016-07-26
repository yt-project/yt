# -*- coding: UTF-8 -*-
#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import contextlib
import mock
import os
import sys
import unittest
import yt.utilities.configure
import yt.config
from yt.config import \
    CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE
from yt.extern.six import StringIO
from yt.extern.six.moves.configparser import NoOptionError, SafeConfigParser


_DUMMY_CFG = ['[yt]', 'loglevel = 49']


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
                yt.utilities.configure.CONFIG = SafeConfigParser()
                if not yt.utilities.configure.CONFIG.has_section('yt'):
                    yt.utilities.configure.CONFIG.add_section('yt')


def tearDownModule():
    for cfgfile in (CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE): 
        if os.path.exists(cfgfile + '.bak_test'):
            os.rename(cfgfile + '.bak_test', cfgfile)


class TestYTConfig(unittest.TestCase):
    def _runYTConfig(self, args):
        args = ['yt-config'] + args
        retcode = 0

        with mock.patch.object(sys, 'argv', args),\
                mock.patch('sys.exit', side_effect=SysExitException) as exit,\
                captureOutput() as output:
            try:
                yt.utilities.configure.main()
            except SysExitException:
                args = exit.mock_calls[0][1]
                retcode = args[0] if len(args) else 0
        return {
            'rc': retcode,
            'stdout': output[0],
            'stderr': output[1]
        }

class TestYTConfigCommands(TestYTConfig):
    def testConfigCommands(self):
        self.assertFalse(os.path.exists(CURRENT_CONFIG_FILE))

        info = self._runYTConfig(['--help'])
        self.assertEqual(info['rc'], 0)
        self.assertEqual(info['stderr'], '')
        self.assertIn('Get and set configuration values for yt',
                      info['stdout'])

        info = self._runYTConfig(['list'])
        self.assertEqual(info['rc'], 0)
        self.assertIn('[yt]', info['stdout'])

        info = self._runYTConfig(['set', 'yt', '__parallel', 'True'])
        self.assertEqual(info['rc'], 0)

        info = self._runYTConfig(['get', 'yt', '__parallel'])
        self.assertEqual(info['rc'], 0)
        self.assertEqual(info['stdout'].strip(), 'True')

        info = self._runYTConfig(['rm', 'yt', '__parallel'])
        self.assertEqual(info['rc'], 0)

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

    def tearDown(self):
        if os.path.exists(CURRENT_CONFIG_FILE):
            os.remove(CURRENT_CONFIG_FILE)
        if os.path.exists(_OLD_CONFIG_FILE + '.bak'):
            os.remove(_OLD_CONFIG_FILE + '.bak')

    def testConfigMigration(self):
        self.assertFalse(os.path.exists(CURRENT_CONFIG_FILE))
        self.assertTrue(os.path.exists(_OLD_CONFIG_FILE))
        
        info = self._runYTConfig(['migrate'])
        self.assertEqual(info['rc'], 0)

        self.assertTrue(os.path.exists(CURRENT_CONFIG_FILE))
        self.assertFalse(os.path.exists(_OLD_CONFIG_FILE))
        self.assertTrue(os.path.exists(_OLD_CONFIG_FILE + '.bak'))

        with open(CURRENT_CONFIG_FILE, 'r') as fh:
            new_cfg = ''.join(fh.readlines())
        self.assertEqual(new_cfg.strip().split('\n'), _DUMMY_CFG)
