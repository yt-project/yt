import argparse
import configparser
import os
import sys

from yt.config import LOCAL_CONFIG_FILE  # noqa: F401
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


def write_config(fd):
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
    write_config()
    print(f"Backing up the old config file: {OLD_CONFIG_FILE}.bak")
    os.rename(OLD_CONFIG_FILE, OLD_CONFIG_FILE + ".bak")


def rm_config(section, option):
    option_path = option.split(".")
    CONFIG.remove(section, *option_path)
    write_config()


def main():
    global LOCAL_CONFIG_FILE, GLOBAL_CONFIG_FILE
    parser = argparse.ArgumentParser(
        description="Get and set configuration values for yt"
    )
    subparsers = parser.add_subparsers(help="sub-command help", dest="cmd")

    get_parser = subparsers.add_parser("get", help="get a config value")
    set_parser = subparsers.add_parser("set", help="set a config value")
    rm_parser = subparsers.add_parser("rm", help="remove a config option")
    subparsers.add_parser("migrate", help="migrate old config file")
    list_parser = subparsers.add_parser("list", help="show all config values")
    print_path_parser = subparsers.add_parser(
        "print-path", help="print the path to the config file"
    )

    get_parser.add_argument("section", help="The section containing the option.")
    get_parser.add_argument("option", help="The option to retrieve.")

    set_parser.add_argument("section", help="The section containing the option.")
    set_parser.add_argument("option", help="The option to set.")
    set_parser.add_argument("value", help="The value to set the option to.")

    rm_parser.add_argument(
        "section", help="The section containing the option to remove."
    )
    rm_parser.add_argument("option", help="The option to remove.")

    for p in (set_parser, get_parser, rm_parser, list_parser, print_path_parser):
        p.add_argument(
            "--local",
            action="store_true",
            help="Use a local configuration file instead of the global one.",
        )

    args = parser.parse_args()

    if "local" in args and args.local:
        if LOCAL_CONFIG_FILE is None:
            LOCAL_CONFIG_FILE = os.path.abspath(os.path.join(".", "yt.toml"))
            with open(LOCAL_CONFIG_FILE, "w") as f:
                f.write("[yt]")
        config_file = LOCAL_CONFIG_FILE
    else:
        config_file = GLOBAL_CONFIG_FILE

    CONFIG.read(config_file)

    if args.cmd == "get":
        print(get_config(args.section, args.option))
    elif args.cmd == "set":
        set_config(args.section, args.option, args.value, config_file)
    elif args.cmd == "list":
        write_config(sys.stdout)
    elif args.cmd == "migrate":
        migrate_config()
    elif args.cmd == "rm":
        rm_config(args.section, args.option, config_file)
    elif args.cmd == "print-path":
        print(config_file)


if __name__ == "__main__":
    main()  # pragma: no cover
