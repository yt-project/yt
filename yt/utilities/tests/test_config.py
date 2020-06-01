import contextlib
import os
import sys
import unittest
import unittest.mock as mock
from io import StringIO

import yt.config
import yt.utilities.command_line
from yt.config import OLD_CONFIG_FILE, YTConfig

GLOBAL_CONFIG_FILE = YTConfig.get_global_config_file()
LOCAL_CONFIG_FILE = YTConfig.get_local_config_file()

_TEST_PLUGIN = "_test_plugin.py"
# NOTE: the normalization of the crazy camel-case will be checked
_DUMMY_CFG_INI = f"""[yt]
logLevel = 49
pluginfilename = {_TEST_PLUGIN}
boolean_stuff = True
chunk_size = 3
"""

_DUMMY_CFG_TOML = f"""[yt]
log_level = 49
plugin_filename = "{_TEST_PLUGIN}"
boolean_stuff = true
chunk_size = 3
"""


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
    for cfgfile in (GLOBAL_CONFIG_FILE, OLD_CONFIG_FILE, LOCAL_CONFIG_FILE):
        if os.path.exists(cfgfile):
            os.rename(cfgfile, cfgfile + ".bak_test")

            if cfgfile == GLOBAL_CONFIG_FILE:
                yt.utilities.configure.CONFIG = YTConfig()
                if not yt.utilities.configure.CONFIG.has_section("yt"):
                    yt.utilities.configure.CONFIG.add_section("yt")


def tearDownModule():
    for cfgfile in (GLOBAL_CONFIG_FILE, OLD_CONFIG_FILE, LOCAL_CONFIG_FILE):
        if os.path.exists(cfgfile + ".bak_test"):
            os.rename(cfgfile + ".bak_test", cfgfile)


class TestYTConfig(unittest.TestCase):
    def _runYTConfig(self, args):
        args = ["yt", "config"] + args
        retcode = 0

        with mock.patch.object(sys, "argv", args), mock.patch(
            "sys.exit", side_effect=SysExitException
        ) as exit, captureOutput() as output:
            try:
                yt.utilities.command_line.run_main()
            except SysExitException:
                args = exit.mock_calls[0][1]
                retcode = args[0] if len(args) else 0
        return {"rc": retcode, "stdout": output[0], "stderr": output[1]}

    def _testKeyValue(self, key, val_set, val_get):
        info = self._runYTConfig(["set", "yt", key, str(val_set)])
        self.assertEqual(info["rc"], 0)

        info = self._runYTConfig(["get", "yt", key])
        self.assertEqual(info["rc"], 0)
        self.assertEqual(info["stdout"].strip(), str(val_get))

        info = self._runYTConfig(["rm", "yt", key])
        self.assertEqual(info["rc"], 0)

    def _testKeyTypeError(self, key, val1, val2, expect_error):
        info = self._runYTConfig(["set", "yt", key, str(val1)])
        self.assertEqual(info["rc"], 0)

        if expect_error:
            with self.assertRaises(TypeError):
                info = self._runYTConfig(["set", "yt", key, str(val2)])
        else:
            info = self._runYTConfig(["set", "yt", key, str(val2)])

        info = self._runYTConfig(["rm", "yt", key])
        self.assertEqual(info["rc"], 0)


class TestYTConfigCommands(TestYTConfig):
    def testConfigCommands(self):
        def remove_spaces_and_breaks(s):
            return "".join(s.split())

        self.assertFalse(os.path.exists(GLOBAL_CONFIG_FILE))

        info = self._runYTConfig(["--help"])
        self.assertEqual(info["rc"], 0)
        self.assertEqual(info["stderr"], "")
        self.assertIn(
            remove_spaces_and_breaks("Get and set configuration values for yt"),
            remove_spaces_and_breaks(info["stdout"]),
        )

        info = self._runYTConfig(["list"])
        self.assertEqual(info["rc"], 0)
        self.assertEqual(info["stdout"], "")

        self._testKeyValue("internals.parallel", True, True)
        self._testKeyValue(
            "test_data_dir", "~/yt-data", os.path.expanduser("~/yt-data")
        )
        self._testKeyValue(
            "test_data_dir", "$HOME/yt-data", os.path.expandvars("$HOME/yt-data")
        )

        with self.assertRaises(KeyError):
            self._runYTConfig(["get", "yt", "foo"])

        # Check TypeErrors are raised when changing the type of an entry
        self._testKeyTypeError("foo.bar", "test", 10, expect_error=True)
        self._testKeyTypeError("foo.bar", "test", False, expect_error=True)

        # Check no type error are raised when *not* changing the type
        self._testKeyTypeError("foo.bar", 10, 20, expect_error=False)
        self._testKeyTypeError("foo.bar", "foo", "bar", expect_error=False)

    def tearDown(self):
        if os.path.exists(GLOBAL_CONFIG_FILE):
            os.remove(GLOBAL_CONFIG_FILE)


class TestYTConfigGlobalLocal(TestYTConfig):
    def setUp(self):
        with open(YTConfig.get_local_config_file(), mode="w") as f:
            f.writelines("[yt]\n")
        with open(YTConfig.get_global_config_file(), mode="w") as f:
            f.writelines("[yt]\n")

    def testAmbiguousConfig(self):
        info = self._runYTConfig(["list"])
        self.assertFalse(len(info["rc"]) == 0)

        for cmd in (["list", "--local"], ["list", "--global"]):
            info = self._runYTConfig(cmd)
            self.assertEqual(info["rc"], 0)


class TestYTConfigMigration(TestYTConfig):
    def setUp(self):
        # TODO: rewrite this with pytest fixtures monkeypatching and tmp_path
        if not os.path.exists(os.path.dirname(OLD_CONFIG_FILE)):
            os.makedirs(os.path.dirname(OLD_CONFIG_FILE))

        with open(OLD_CONFIG_FILE, "w") as fh:
            fh.write(_DUMMY_CFG_INI)

        if os.path.exists(GLOBAL_CONFIG_FILE):
            os.remove(GLOBAL_CONFIG_FILE)

    def tearDown(self):
        # TODO: rewrite this with pytest fixtures monkeypatching and tmp_path
        if os.path.exists(GLOBAL_CONFIG_FILE):
            os.remove(GLOBAL_CONFIG_FILE)
        if os.path.exists(OLD_CONFIG_FILE + ".bak"):
            os.remove(OLD_CONFIG_FILE + ".bak")

    def testConfigMigration(self):
        self.assertFalse(os.path.exists(GLOBAL_CONFIG_FILE))
        self.assertTrue(os.path.exists(OLD_CONFIG_FILE))

        info = self._runYTConfig(["migrate"])
        self.assertEqual(info["rc"], 0)

        self.assertTrue(os.path.exists(GLOBAL_CONFIG_FILE))
        self.assertFalse(os.path.exists(OLD_CONFIG_FILE))
        self.assertTrue(os.path.exists(OLD_CONFIG_FILE + ".bak"))

        with open(GLOBAL_CONFIG_FILE, mode="r") as fh:
            new_cfg = fh.read()

        self.assertEqual(new_cfg, _DUMMY_CFG_TOML)
