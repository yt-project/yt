# -*- coding: UTF-8 -*-
#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import shutil
import sys
import argparse
from yt.config import CURRENT_CONFIG_FILE, _OLD_CONFIG_FILE, YTConfigParser

CONFIG = YTConfigParser()
CONFIG.read([CURRENT_CONFIG_FILE])


def get_config(section, option):
    return CONFIG.get(section, option)


def set_config(section, option, value):
    if not CONFIG.has_section(section):
        CONFIG.add_section(section)
    CONFIG.set(section, option, value)
    write_config()


def write_config(fd=None):
    if fd is None:
        with open(CURRENT_CONFIG_FILE, 'w') as fd:
            CONFIG.write(fd)
    else:
        CONFIG.write(fd)

def migrate_config():
    if not os.path.exists(_OLD_CONFIG_FILE):
        print('Old config not found.')
        sys.exit()
    CONFIG.read(_OLD_CONFIG_FILE)
    print('Writing a new config file to: {}'.format(CURRENT_CONFIG_FILE))
    write_config()
    print('Backing up the old config file: {}.bak'.format(_OLD_CONFIG_FILE))
    os.rename(_OLD_CONFIG_FILE, _OLD_CONFIG_FILE + '.bak')

    old_config_dir = os.path.dirname(_OLD_CONFIG_FILE)
    plugin_file = CONFIG.get('yt', 'pluginfilename')
    if plugin_file and \
            os.path.exists(os.path.join(old_config_dir, plugin_file)):
        print('Migrating plugin file {} to new location'.format(plugin_file))
        shutil.copyfile(
            os.path.join(old_config_dir, plugin_file),
            os.path.join(os.path.dirname(CURRENT_CONFIG_FILE), plugin_file))
        print('Backing up the old plugin file: {}.bak'.format(_OLD_CONFIG_FILE))
        plugin_file = os.path.join(old_config_dir, plugin_file)
        os.rename(plugin_file, plugin_file + '.bak')


def rm_config(section, option):
    CONFIG.remove_option(section, option)
    write_config()


def main():
    parser = argparse.ArgumentParser(
        description='Get and set configuration values for yt')
    subparsers = parser.add_subparsers(help='sub-command help', dest='cmd')

    get_parser = subparsers.add_parser('get', help='get a config value')
    set_parser = subparsers.add_parser('set', help='set a config value')
    rm_parser = subparsers.add_parser('rm', help='remove a config option')
    subparsers.add_parser('migrate', help='migrate old config file')
    subparsers.add_parser('list', help='show all config values')

    get_parser.add_argument(
        'section', help='The section containing the option.')
    get_parser.add_argument('option', help='The option to retrieve.')

    set_parser.add_argument(
        'section', help='The section containing the option.')
    set_parser.add_argument('option', help='The option to set.')
    set_parser.add_argument('value', help='The value to set the option to.')

    rm_parser.add_argument(
        'section', help='The section containing the option to remove.')
    rm_parser.add_argument('option', help='The option to remove.')

    args = parser.parse_args()

    if args.cmd == 'get':
        print(get_config(args.section, args.option))
    elif args.cmd == 'set':
        set_config(args.section, args.option, args.value)
    elif args.cmd == 'list':
        write_config(sys.stdout)
    elif args.cmd == 'migrate':
        migrate_config()
    elif args.cmd == 'rm':
        rm_config(args.section, args.option)

if __name__ == '__main__':
    main()  # pragma: no cover
