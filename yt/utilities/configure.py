import configparser
import os
import sys

from yt.config import GLOBAL_CONFIG_FILE, OLD_CONFIG_FILE, YTConfig, ytcfg_defaults

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


def migrate_config():
    if not os.path.exists(OLD_CONFIG_FILE):
        print("Old config not found.")
        sys.exit(1)

    old_config = configparser.RawConfigParser()
    # Preserve case:
    # See https://stackoverflow.com/questions/1611799/preserve-case-in-configparser
    old_config.optionxform = str
    old_config.read(OLD_CONFIG_FILE)

    default_keys = {k.lower(): k for k in ytcfg_defaults["yt"].keys()}

    config_as_dict = {}
    for section in old_config:
        if section == "DEFAULT":
            continue
        config_as_dict[section] = {}
        for key, value in old_config[section].items():
            # Cast value to the most specific type possible
            cast_value = _cast_value_helper(value)

            # Normalize the key (if present in the defaults)
            if key.lower() in default_keys and section == "yt":
                normalized_key = default_keys[key.lower()]
            else:
                normalized_key = key

            config_as_dict[section][normalized_key] = cast_value

    CONFIG.update(config_as_dict)

    print(f"Writing a new config file to: {GLOBAL_CONFIG_FILE}")
    write_config(GLOBAL_CONFIG_FILE)
    print(f"Backing up the old config file: {OLD_CONFIG_FILE}.bak")
    os.rename(OLD_CONFIG_FILE, OLD_CONFIG_FILE + ".bak")


def rm_config(section, option, config_file):
    option_path = option.split(".")
    CONFIG.remove(section, *option_path)
    write_config(config_file)
