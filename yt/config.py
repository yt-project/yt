import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import toml

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
CURRENT_CONFIG_FILE = os.path.join(CONFIG_DIR, "yt.toml")

if os.path.exists(OLD_CONFIG_FILE):
    if os.path.exists(CURRENT_CONFIG_FILE):
        msg = (
            "The configuration file {} is deprecated. "
            "Please manually remove it to suppress this warning."
        )
        warnings.warn(msg.format(OLD_CONFIG_FILE, CURRENT_CONFIG_FILE))
    else:
        msg = (
            "The configuration file {} is deprecated. "
            "Please migrate your config to {} by running: "
            "'yt config migrate'"
        )
        warnings.warn(msg.format(OLD_CONFIG_FILE, CURRENT_CONFIG_FILE))
        sys.exit(1)


if not os.path.exists(CURRENT_CONFIG_FILE):
    cfg = {"yt": {}}
    try:
        with open(CURRENT_CONFIG_FILE, mode="w") as fd:
            toml.dump(cfg, fd)
    except IOError:
        warnings.warn("unable to write new config file")


class YTConfig:
    def __init__(self, defaults=None):
        self.defaults = defaultdict(dict)
        self.defaults.update(defaults)

        self.values = defaultdict(dict)
        self.update(defaults)

    @staticmethod
    def _get_helper(
        root, option, strict=True, lower_case=False, fallback=None, use_fallback=False
    ):
        *option_path, option_name = option
        # This works as follow: if we try to access
        # field > gas > density > lognorm
        # we try in this order:
        #   field > gas > density > lognorm
        #   field > gas > lognorm
        #   field > lognorm

        node = {}
        first_pass = True
        while len(option_path) > 0 or first_pass:
            first_pass = False
            try:
                node = root
                for k in option_path:
                    if lower_case:
                        lower_keys = {k.lower(): k for k in node.keys()}
                        k_to_load = lower_keys[k.lower()]
                    else:
                        k_to_load = k
                    node = node[k_to_load]

                if lower_case:
                    lower_keys = {k.lower(): k for k in node.keys()}
                    option_name_lower = lower_keys[option_name.lower()]
                    return option_name_lower, node[option_name_lower]
                else:
                    return option_name, node[option_name]
            except KeyError as e:
                if strict and not use_fallback:
                    raise e
                else:
                    # Shorten the path and try again
                    option_path = option_path[:-1]

        if use_fallback:
            return fallback
        else:
            option_as_str = ".".join(option_path + [option_name])
            raise KeyError(f"Could not find {option_as_str} in configuration.")

    def get(self, section, *option, strict=True, **kwargs):
        _, value = self._get_helper(
            self.values, [section] + list(option), strict=strict, **kwargs
        )
        return value

    def update(self, new_values):
        def copy_helper(dict_a, dict_b):
            # Copies entries form dict_a in dict_b, inplace
            for k in dict_a.keys():
                if k not in dict_b:
                    dict_b[k] = dict_a[k]
                    continue
                if isinstance(dict_a[k], dict):
                    copy_helper(dict_a[k], dict_b[k])
                else:
                    dict_b[k] = dict_a[k]

        copy_helper(new_values, self.values)

    def has_section(self, section):
        return section in self.values

    def add_section(self, section):
        if not self.has_section(section):
            self.values[section] = {}

    def remove_section(self, section):
        if self.has_section(section):
            self.values.pop(section)

    def __setitem__(self, key, value):
        *option_path, option_name = key

        _, node = self._get_helper(self.values, option_path, strict=True)

        # Get the corresponding defaults
        default_option_name, default_value = self._get_helper(
            self.defaults,
            key,
            strict=True,
            fallback=None,
            use_fallback=True,
            lower_case=True,
        )

        # ... and check the type match
        if default_value is not None:
            default_type = type(default_value)
            if not type(value) == default_type:
                raise TypeError(
                    (
                        "Error when setting %s to %s. "
                        "Expected a value with type %s, got instead a value of type %s."
                    )
                    % (key, value, default_type, type(value))
                )

        # ... and the case match
        if default_option_name != option_name:
            raise KeyError(
                "Error when setting %s to %s. Did you mean `%s`?"
                % (key, value, list(key[:-1]) + [default_option_name])
            )
        node[option_name] = value

    def set(self, section, option, value):
        if not isinstance(option, (tuple, list)):
            option = (option,)
        self[(section, *option)] = value

    def read(self, file_names):
        if not isinstance(file_names, (tuple, list)):
            file_names = (file_names,)

        file_names_read = []
        for fname in file_names:
            if not os.path.exists(fname):
                continue
            self.update(toml.load(fname))
            file_names_read.append(fname)

        return file_names_read

    def write(self, fd):
        def cleaner(d):
            new_d = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    if len(v) > 0:
                        new_d[k] = cleaner(v)
                else:
                    new_d[k] = v
            return new_d

        # Clean up the sections
        cleaned_values = cleaner(self.values)
        config_as_str = toml.dumps(cleaned_values)

        try:
            fd.write(config_as_str)
        except AttributeError:
            with open(fd, mode="w") as fdd:
                fdd.write(config_as_str)


# Walk the tree up until we find a config file
ytcfg = YTConfig(defaults=ytcfg_defaults)
if os.path.exists(CURRENT_CONFIG_FILE):
    ytcfg.read(CURRENT_CONFIG_FILE)
cwd = Path(".").absolute()
while cwd.parent != cwd:
    cfg_file = cwd / "yt.toml"
    if cfg_file.exists():
        ytcfg.read(cfg_file)
        break
    cwd = cwd.parent
