import os
import shutil
import tempfile
import unittest

import yt
from yt.config import ytcfg
from yt.testing import assert_raises, fake_random_ds
from yt.utilities.configure import YTConfig, config_dir

_TEST_PLUGIN = "_test_plugin.py"

_DUMMY_CFG_TOML = f"""[yt]
log_level = 49
plugin_filename = "{_TEST_PLUGIN}"
boolean_stuff = true
chunk_size = 3
"""

TEST_PLUGIN_FILE = """
import numpy as np

def _myfunc(field, data):
    return np.random.random(data[('gas', 'density')].shape)
add_field(('gas', 'random'), dimensions='dimensionless',
          function=_myfunc, units='auto', sampling_type='local')
constant = 3
def myfunc():
    return constant*4
foobar = 17
"""


class TestPluginFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
        cls.tmpdir = tempfile.mkdtemp()
        os.mkdir(os.path.join(cls.tmpdir, "yt"))
        os.environ["XDG_CONFIG_HOME"] = cls.tmpdir
        with open(YTConfig.get_global_config_file(), mode="w") as fh:
            fh.write(_DUMMY_CFG_TOML)
        cls.plugin_path = os.path.join(config_dir(), ytcfg.get("yt", "plugin_filename"))
        with open(cls.plugin_path, mode="w") as fh:
            fh.write(TEST_PLUGIN_FILE)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        if cls.xdg_config_home:
            os.environ["XDG_CONFIG_HOME"] = cls.xdg_config_home
        else:
            os.environ.pop("XDG_CONFIG_HOME")

    def testCustomField(self):
        msg = f"INFO:yt:Loading plugins from {self.plugin_path}"
        with self.assertLogs("yt", level="INFO") as cm:
            yt.enable_plugins()
            self.assertEqual(cm.output, [msg])

        ds = fake_random_ds(16)
        dd = ds.all_data()
        self.assertEqual(str(dd[("gas", "random")].units), "dimensionless")
        self.assertEqual(dd[("gas", "random")].shape, dd[("gas", "density")].shape)
        assert yt.myfunc() == 12
        assert_raises(AttributeError, getattr, yt, "foobar")
