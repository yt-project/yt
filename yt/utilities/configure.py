import os

from yt.config import YTConfig

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
