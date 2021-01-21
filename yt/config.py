import os
import warnings

import toml
from more_itertools import always_iterable

from yt.utilities.configuration_tree import ConfigNode

ytcfg_defaults = {}

ytcfg_defaults["yt"] = dict(
    serialize=False,
    only_deserialize=False,
    time_functions=False,
    log_file=False,
    colored_logs=False,
    suppress_stream_logging=False,
    stdout_stream_logging=False,
    log_level=20,
    inline=False,
    num_threads=-1,
    store_parameter_files=False,
    parameter_file_store="parameter_files.csv",
    maximum_stored_datasets=500,
    skip_dataset_cache=True,
    load_field_plugins=False,
    plugin_filename="my_plugins.py",
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
        within_testing=False,
        within_pytest=False,
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


class YTConfig:
    def __init__(self, defaults=None):
        if defaults is None:
            defaults = {}
        self.config_root = ConfigNode(None)

    def get(self, section, *keys, callback=None, **kwargs):
        if callback is None:

            def callback(leaf):
                return leaf.value

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
            [section] + list(keys), value, extra_data=metadata
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
        file_names_read = []
        for fname in always_iterable(file_names):
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

    @staticmethod
    def get_global_config_file():
        return os.path.join(CONFIG_DIR, "yt.toml")

    @staticmethod
    def get_local_config_file():
        return os.path.join(os.path.abspath(os.curdir), "yt.toml")


OLD_CONFIG_FILE = os.path.join(CONFIG_DIR, "ytrc")
_global_config_file = YTConfig.get_global_config_file()
_local_config_file = YTConfig.get_local_config_file()

if os.path.exists(OLD_CONFIG_FILE):
    if os.path.exists(_global_config_file):
        msg = (
            f"The configuration file {OLD_CONFIG_FILE} is deprecated in "
            f"favor of {_global_config_file}. Currently, both are present. "
            "Please manually remove the deprecated one to silence "
            "this warning."
        )
        warnings.warn(msg)
    else:
        # We have an issue here: when calling from the command line,
        # we do not want this to exit, as it would prevent `yt config migrate`
        # from running. The issue is that yt.config (this file) is imported
        # from yt.__init__, which is imported *before* yt.utilities.configure
        # is executed, so the latter cannot set some internal variable before
        # we arrive here.
        # The workaround here relies on inspecting the call stack and hopefully
        # detect we were called from the CLI.
        import inspect

        stack = inspect.stack()
        if len(stack) < 2 or stack[-2].function != "importlib_load_entry_point":
            msg = (
                f"The configuration file {OLD_CONFIG_FILE} is deprecated. "
                f"Please migrate your config to {_global_config_file} by running: "
                "'yt config migrate'"
            )
            raise SystemExit(msg)


if not os.path.exists(_global_config_file):
    cfg = {"yt": {}}
    try:
        with open(_global_config_file, mode="w") as fd:
            toml.dump(cfg, fd)
    except OSError:
        warnings.warn("unable to write new config file")


# Load the config
ytcfg = YTConfig()
ytcfg.update(ytcfg_defaults, metadata={"source": "defaults"})

# Try loading the local config first, otherwise fall back to global config
if os.path.exists(_local_config_file):
    ytcfg.read(_local_config_file)
elif os.path.exists(_global_config_file):
    ytcfg.read(_global_config_file)
