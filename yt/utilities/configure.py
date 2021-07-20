import configparser
import os
import sys

from yt.config import YTConfig, old_config_file, ytcfg_defaults

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
    if not os.path.exists(old_config_file()):
        print("Old config not found.")
        sys.exit(1)

    old_config = configparser.RawConfigParser()
    # Preserve case:
    # See https://stackoverflow.com/questions/1611799/preserve-case-in-configparser
    old_config.optionxform = str
    old_config.read(old_config_file())

    # In order to migrate, we'll convert everything to lowercase, and map that
    # to the new snake_case convention
    def normalize_key(key):
        return key.replace("_", "").lower()

    def usesCamelCase(key):
        if key != key.lower():
            return True
        else:
            return False

    old_keys_to_new = {normalize_key(k): k for k in ytcfg_defaults["yt"].keys()}

    config_as_dict = {}
    for section in old_config:
        if section == "DEFAULT":
            continue
        config_as_dict[section] = {}
        for key, value in old_config[section].items():
            # Cast value to the most specific type possible
            cast_value = _cast_value_helper(value)

            # Normalize the key (if present in the defaults)
            if normalize_key(key) in old_keys_to_new and section == "yt":
                new_key = old_keys_to_new[normalize_key(key)]
            else:
                new_key = key

            config_as_dict[section][new_key] = cast_value

    CONFIG.update(config_as_dict)

    global_config_file = YTConfig.get_global_config_file()
    print(f"Writing a new config file to: {global_config_file}")
    write_config(global_config_file)
    print(f"Backing up the old config file: {old_config_file()}.bak")
    os.rename(old_config_file(), old_config_file() + ".bak")


def rm_config(section, option, config_file):
    option_path = option.split(".")
    CONFIG.remove(section, *option_path)
    write_config(config_file)
