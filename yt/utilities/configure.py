import argparse
import configparser
import os
import sys

from yt.config import CURRENT_CONFIG_FILE, OLD_CONFIG_FILE, YTConfig

CONFIG = YTConfig()
CONFIG.read(CURRENT_CONFIG_FILE)


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


def set_config(section, option, value):
    if not CONFIG.has_section(section):
        CONFIG.add_section(section)

    option_path = option.split(".")
    CONFIG.set(section, option_path, _cast_value_helper(value))
    write_config()


def write_config(fd=None):
    if fd is None:
        CONFIG.write(CURRENT_CONFIG_FILE)
    else:
        CONFIG.write(fd)


def migrate_config():
    if not os.path.exists(OLD_CONFIG_FILE):
        print("Old config not found.")
        sys.exit(1)

    old_config = configparser.RawConfigParser()
    # Preserve case:
    # See https://stackoverflow.com/questions/1611799/preserve-case-in-configparser
    old_config.optionxform = str
    old_config.read(OLD_CONFIG_FILE)

    config_as_dict = {}
    for section in old_config:
        if section == "DEFAULT":
            continue
        config_as_dict[section] = {}
        for key, value in old_config[section].items():
            config_as_dict[section][key] = _cast_value_helper(value)

    print(config_as_dict)

    CONFIG.update(config_as_dict)

    print(f"Writing a new config file to: {CURRENT_CONFIG_FILE}")
    write_config()
    print(f"Backing up the old config file: {OLD_CONFIG_FILE}.bak")
    os.rename(OLD_CONFIG_FILE, OLD_CONFIG_FILE + ".bak")


def rm_config(section, option):
    *option_path, option_name = option.split(".")
    parent = CONFIG.get(section, *option_path)
    parent.pop(option_name)
    write_config()


def main():
    parser = argparse.ArgumentParser(
        description="Get and set configuration values for yt"
    )
    subparsers = parser.add_subparsers(help="sub-command help", dest="cmd")

    get_parser = subparsers.add_parser("get", help="get a config value")
    set_parser = subparsers.add_parser("set", help="set a config value")
    rm_parser = subparsers.add_parser("rm", help="remove a config option")
    subparsers.add_parser("migrate", help="migrate old config file")
    subparsers.add_parser("list", help="show all config values")

    get_parser.add_argument("section", help="The section containing the option.")
    get_parser.add_argument("option", help="The option to retrieve.")

    set_parser.add_argument("section", help="The section containing the option.")
    set_parser.add_argument("option", help="The option to set.")
    set_parser.add_argument("value", help="The value to set the option to.")

    rm_parser.add_argument(
        "section", help="The section containing the option to remove."
    )
    rm_parser.add_argument("option", help="The option to remove.")

    args = parser.parse_args()

    if args.cmd == "get":
        print(get_config(args.section, args.option))
    elif args.cmd == "set":
        set_config(args.section, args.option, args.value)
    elif args.cmd == "list":
        write_config(sys.stdout)
    elif args.cmd == "migrate":
        migrate_config()
    elif args.cmd == "rm":
        rm_config(args.section, args.option)


if __name__ == "__main__":
    main()  # pragma: no cover
