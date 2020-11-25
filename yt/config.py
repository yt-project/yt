import os
import sys
import warnings
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
        if defaults is None:
            defaults = {"yt": {}}
        self.defaults = defaults
        self.values = {}
        self.update(defaults)

    def get(self, section, *path, strict=True, **kwargs):
        config = self.values[section]

        # This works as follow: if we try to access
        # field > gas > density > lognorm
        # we try in this order:
        #   field > gas > density > lognorm
        #   field > gas > lognorm
        #   field > lognorm
        if len(path) == 0:
            return config

        ok = False
        node = None
        use_fallback = "fallback" in kwargs
        fallback = kwargs.pop("fallback", None)
        while len(path) > 0:
            try:
                node = config
                for k in path:
                    node = node[k]
                ok = True
                break
            except KeyError as e:
                if strict and not use_fallback:
                    raise e
                else:
                    path = path[:-1]

        if not ok and use_fallback:
            return fallback
        elif not ok:
            raise KeyError(f"Could not find {section}, {path} in configuration.")

        return node

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

    def __setitem__(self, key, value):
        section, *option_path, option_name = key
        if not self.has_section(section):
            raise KeyError

        node = self.values[section]
        for p in option_path:
            if p not in node:
                node[p] = {}
            node = node[p]

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
ytcfg = YTConfig(ytcfg_defaults)
if os.path.exists(CURRENT_CONFIG_FILE):
    ytcfg.read(CURRENT_CONFIG_FILE)
cwd = Path(".").absolute()
while cwd.parent != cwd:
    cfg_file = cwd / "yt.toml"
    if cfg_file.exists():
        ytcfg.read(cfg_file)
        break
    cwd = cwd.parent

# Now we have parsed the config file.  Overrides come from the command line.

# This should be implemented at some point.  The idea would be to have a set of
# command line options, fed through an option parser, that would override
# the settings in ytcfg.  *However*, because we want to have the command-line
# scripts work, we'd probably want to have them only be long options, and also
# along the lines of --yt-something=somethingelse.  The command line scripts
# would then not get their options from sys.argv, but instead from this module.
