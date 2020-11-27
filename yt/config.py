import os
import sys
import warnings
from itertools import chain
from pathlib import Path

import toml

from yt.utilities.configuration_tree import ConfigNode

ytcfg_defaults = {}

ytcfg_defaults["yt"] = dict(
    serialize=False,
    only_deserialize=False,
    timefunctions=False,
    logFile=False,
    coloredLogs=False,
    suppressStreamLogging=False,
    stdoutStreamLogging=False,
    loglevel=20,
    inline=False,
    numthreads=-1,
    storeParameterFiles=False,
    parameterFileStore="parameter_files.csv",
    maximumStoredDatasets=500,
    skip_dataset_cache=True,
    loadFieldPlugins=False,
    pluginFilename="my_plugins.py",
    parallel_traceback=False,
    pasteboard_repo="",
    reconstruct_index=True,
    test_storage_dir="/does/not/exist",
    test_data_dir="/does/not/exist",
    enzo_db="",
    hub_url="https://girder.hub.yt/api/v1",
    hub_api_key="",
    hub_sandbox="/collection/yt_sandbox/data",
    notebook_password="",
    answer_testing_tolerance=3,
    answer_testing_bitwise=False,
    gold_standard_filename="gold311",
    local_standard_filename="local001",
    answer_tests_url="http://answers.yt-project.org/{1}_{2}",
    sketchfab_api_key="None",
    imagebin_api_key="e1977d9195fe39e",
    imagebin_upload_url="https://api.imgur.com/3/upload",
    imagebin_delete_url="https://api.imgur.com/3/image/{delete_hash}",
    curldrop_upload_url="http://use.yt/upload",
    thread_field_detection=False,
    ignore_invalid_unit_operation_errors=False,
    chunk_size=1000,
    xray_data_dir="/does/not/exist",
    supp_data_dir="/does/not/exist",
    default_colormap="arbre",
    ray_tracing_engine="embree",
    internals=dict(
        withintesting=False,
        withinpytest=False,
        parallel=False,
        strict_requires=False,
        global_parallel_rank=0,
        global_parallel_size=1,
        topcomm_parallel_rank=0,
        topcomm_parallel_size=1,
        command_line=False,
    ),
)

CONFIG_DIR = os.environ.get(
    "XDG_CONFIG_HOME", os.path.join(os.path.expanduser("~"), ".config", "yt")
)
if not os.path.exists(CONFIG_DIR):
    try:
        os.makedirs(CONFIG_DIR)
    except OSError:
        warnings.warn("unable to create yt config directory")

OLD_CONFIG_FILE = os.path.join(CONFIG_DIR, "ytrc")
GLOBAL_CONFIG_FILE = os.path.join(CONFIG_DIR, "yt.toml")
LOCAL_CONFIG_FILE = None

if os.path.exists(OLD_CONFIG_FILE):
    if os.path.exists(GLOBAL_CONFIG_FILE):
        msg = (
            "The configuration file {} is deprecated. "
            "Please manually remove it to suppress this warning."
        )
        warnings.warn(msg.format(OLD_CONFIG_FILE, GLOBAL_CONFIG_FILE))
    else:
        msg = (
            "The configuration file {} is deprecated. "
            "Please migrate your config to {} by running: "
            "'yt config migrate'"
        )
        warnings.warn(msg.format(OLD_CONFIG_FILE, GLOBAL_CONFIG_FILE))
        sys.exit(1)


if not os.path.exists(GLOBAL_CONFIG_FILE):
    cfg = {"yt": {}}
    try:
        with open(GLOBAL_CONFIG_FILE, mode="w") as fd:
            toml.dump(cfg, fd)
    except IOError:
        warnings.warn("unable to write new config file")


class YTConfig:
    def __init__(self, defaults=None):
        if defaults is None:
            defaults = {}
        self.config_root = ConfigNode(None)

    def get(self, section, *keys, callback=None, **kwargs):
        if callback is None:
            callback = lambda leaf: leaf.value  # noqa: E731

        return self.config_root.get_leaf(section, *keys, callback=callback)

    def get_most_specific(self, section, *keys, **kwargs):
        use_fallback = "fallback" in kwargs
        fallback = kwargs.pop("fallback", None)
        try:
            return self.config_root.get_deepest_leaf(section, *keys)
        except KeyError as err:
            if use_fallback:
                return fallback
            else:
                raise err

    def update(self, new_values, metadata=None):
        if metadata is None:
            metadata = {}
        self.config_root.update(new_values, metadata)

    def has_section(self, section):
        try:
            self.config_root.get_child(section)
            return True
        except KeyError:
            return False

    def add_section(self, section):
        self.config_root.add_child(section)

    def remove_section(self, section):
        if self.has_section(section):
            self.config_root.remove_child(section)
            return True
        else:
            return False

    def set(self, *args, metadata=None):
        section, *keys, value = args
        if metadata is None:
            metadata = {"source": "runtime"}
        self.config_root.upsert_from_list(
            [section] + list(keys), value, extraData=metadata
        )

    def __setitem__(self, args, value):
        section, *keys = args
        self.set(section, *keys, value, metadata=None)

    def __getitem__(self, key):
        section, *keys = key
        return self.get(section, *keys)

    def remove(self, *args):
        self.config_root.pop_leaf(args)

    def read(self, file_names):
        if not isinstance(file_names, (tuple, list)):
            file_names = (file_names,)

        file_names_read = []
        for fname in file_names:
            if not os.path.exists(fname):
                continue
            metadata = {"source": f"file: {fname}"}
            self.update(toml.load(fname), metadata=metadata)
            file_names_read.append(fname)

        return file_names_read

    def write(self, file_handler):
        value = self.config_root.as_dict()
        config_as_str = toml.dumps(value)

        try:
            # Assuming file_handler has a write attribute
            file_handler.write(config_as_str)
        except AttributeError:
            # Otherwise we expect a path to a file
            with open(file_handler, mode="w") as fh:
                fh.write(config_as_str)


# Walk the tree up until we find a config file
ytcfg = YTConfig()
ytcfg.update(ytcfg_defaults, metadata={"source": "defaults"})

if os.path.exists(GLOBAL_CONFIG_FILE):
    ytcfg.read(GLOBAL_CONFIG_FILE)

cwd = Path.cwd().resolve()
for folder in chain([cwd], cwd.parents):
    cfg_file = folder / "yt.toml"
    if cfg_file.exists():
        ytcfg.read(cfg_file)
        LOCAL_CONFIG_FILE = str(cfg_file)
        break
