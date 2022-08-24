import os
import sys
import warnings
from pathlib import Path
from typing import Callable, List

import tomli_w
from more_itertools import always_iterable

from yt.utilities.configuration_tree import ConfigLeaf, ConfigNode

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

configuration_callbacks: List[Callable[["YTConfig"], None]] = []


def config_dir():
    config_root = os.environ.get(
        "XDG_CONFIG_HOME", os.path.join(os.path.expanduser("~"), ".config")
    )
    conf_dir = os.path.join(config_root, "yt")
    return conf_dir


class YTConfig:
    def __init__(self, defaults=None):
        if defaults is None:
            defaults = {}
        self.config_root = ConfigNode(None)

    def get(self, section, *keys, callback=None):
        node_or_leaf = self.config_root.get(section, *keys)
        if isinstance(node_or_leaf, ConfigLeaf):
            if callback is not None:
                return callback(node_or_leaf)
            return node_or_leaf.value
        return node_or_leaf

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

    def remove(self, *args):
        self.config_root.pop_leaf(args)

    def read(self, file_names):
        file_names_read = []
        for fname in always_iterable(file_names):
            if not os.path.exists(fname):
                continue
            metadata = {"source": f"file: {fname}"}
            try:
                with open(fname, "rb") as fh:
                    data = tomllib.load(fh)
            except tomllib.TOMLDecodeError as exc:
                warnings.warn(
                    f"Could not load configuration file {fname} (invalid TOML: {exc})"
                )
            else:
                self.update(data, metadata=metadata)
                file_names_read.append(fname)

        return file_names_read

    def write(self, file_handler):
        value = self.config_root.as_dict()
        config_as_str = tomli_w.dumps(value)

        try:
            file_path = Path(file_handler)
        except TypeError:
            if not hasattr(file_handler, "write"):
                raise TypeError(
                    f"Expected a path to a file, or a writable object, got {file_handler}"
                ) from None
            file_handler.write(config_as_str)
        else:
            pdir = file_path.parent
            if not pdir.exists():
                warnings.warn(f"{pdir!s} does not exist, creating it (recursively)")
                os.makedirs(pdir)
            file_path.write_text(config_as_str)

    @staticmethod
    def get_global_config_file():
        return os.path.join(config_dir(), "yt.toml")

    @staticmethod
    def get_local_config_file():
        return os.path.join(os.path.abspath(os.curdir), "yt.toml")

    def __setitem__(self, args, value):
        section, *keys = always_iterable(args)
        self.set(section, *keys, value, metadata=None)

    def __getitem__(self, key):
        section, *keys = always_iterable(key)
        return self.get(section, *keys)

    def __contains__(self, item):
        return item in self.config_root

    # Add support for IPython rich display
    # see https://ipython.readthedocs.io/en/stable/config/integrating.html
    def _repr_json_(self):
        return self.config_root._repr_json_()


CONFIG = YTConfig()


def _cast_bool_helper(value):
    if value == "True":
        return True
    elif value == "False":
        return False
    else:
        raise ValueError("Cannot safely cast to bool")


def _expand_all(s):
    return os.path.expandvars(os.path.expanduser(s))


def _cast_value_helper(value, types=(_cast_bool_helper, int, float, _expand_all)):
    for t in types:
        try:
            retval = t(value)
            return retval
        except ValueError:
            pass


def get_config(section, option):
    *option_path, option_name = option.split(".")
    return CONFIG.get(section, *option_path, option_name)


def set_config(section, option, value, config_file):
    if not CONFIG.has_section(section):
        CONFIG.add_section(section)

    option_path = option.split(".")
    CONFIG.set(section, *option_path, _cast_value_helper(value))
    write_config(config_file)


def write_config(config_file):
    CONFIG.write(config_file)


def rm_config(section, option, config_file):
    option_path = option.split(".")
    CONFIG.remove(section, *option_path)
    write_config(config_file)
