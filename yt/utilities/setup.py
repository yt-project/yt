#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import os.path
import glob
import platform


# snatched from PyTables
def add_from_path(envname, dirs):
    try:
        dirs.extend(os.environ[envname].split(os.pathsep))
    except KeyError:
        pass


# snatched from PyTables
def add_from_flags(envname, flag_key, dirs):
    for flag in os.environ.get(envname, "").split():
        if flag.startswith(flag_key):
            dirs.append(flag[len(flag_key):])


# snatched from PyTables
def get_default_dirs():
    default_header_dirs = []
    default_library_dirs = []

    add_from_path("CPATH", default_header_dirs)
    add_from_path("C_INCLUDE_PATH", default_header_dirs)
    add_from_flags("CPPFLAGS", "-I", default_header_dirs)
    default_header_dirs.extend(
        ['/usr/include', '/usr/local/include', '/usr/X11']
    )

    _archs = ['lib64', 'lib']
    if platform.system() == 'Linux':
        distname, version, did = platform.linux_distribution()
        if distname.lower() in ('ubuntu', 'debian'):
            _archs.extend(
                ['lib/x86_64-linux-gnu',
                 'lib/i686-linux-gnu',
                 'lib/i386-linux-gnu']
            )

    add_from_flags("LDFLAGS", "-L", default_library_dirs)
    default_library_dirs.extend(
        os.path.join(_tree, _arch)
        for _tree in ('/usr', '/usr/local', '/usr/X11', '/')
        for _arch in _archs
    )
    return default_header_dirs, default_library_dirs


def get_location_from_env(env):
    env_dir = os.environ[env]
    env_inc = os.path.join(env_dir, "include")
    env_lib = os.path.join(env_dir, "lib")
    print("%s_LOCATION: %s: %s, %s"
          % (env.split('_')[0], env, env_inc, env_lib))
    return (env_inc, env_lib)


def get_location_from_cfg(cfg):
    cfg_dir = open(cfg).read().strip()
    cfg_inc = os.path.join(cfg_dir, "include")
    cfg_lib = os.path.join(cfg_dir, "lib")
    print("%s_LOCATION: %s: %s, %s"
          % (cfg.split('.')[0].upper(), cfg, cfg_inc, cfg_lib))
    return (cfg_inc, cfg_lib)


def check_prefix(inc_dir, lib_dir):
    if platform.system() == 'Linux':
        distname, version, did = platform.linux_distribution()
        if distname.lower() in ('ubuntu', 'debian'):
            print("Since you are using multiarch distro it's hard to detect")
            print("whether library matches the header file. We will assume")
            print("it does. If you encounter any build failures please use")
            print("proper cfg files to provide path to the dependencies")
            print("")
            return (inc_dir, lib_dir)
    prefix = os.path.commonprefix([inc_dir, lib_dir]).rstrip('/\\')
    if prefix is not '' and prefix == os.path.dirname(inc_dir):
        return (inc_dir, lib_dir)
    else:
        print("It seems that include prefix is different from lib prefix")
        print("Please use either env variable or cfg to set proper path")
        return (None, None)


def get_location_from_ctypes(header, library):
    yt_inst = os.environ.get('YT_DEST')
    if yt_inst is not None:
        # since we prefer installation via script, make sure
        # that YT_DEST path take precedence above all else
        return (os.path.join(yt_inst, 'include'), os.path.join(yt_inst, 'lib'))

    try:
        import ctypes
        import ctypes.util
    except ImportError:
        return (None, None)

    target_inc, target_libdir = None, None
    default_header_dirs, default_library_dirs = get_default_dirs()
    for inc_prefix in default_header_dirs:
        if os.path.isfile(os.path.join(inc_prefix, header)):
            target_inc = inc_prefix

    target_libfile = ctypes.util.find_library(library)
    if None in (target_inc, target_libfile):
        # either header or lib was not found, abort now
        return (None, None)
    if os.path.isfile(target_libfile):
        return check_prefix(target_inc, os.path.dirname(target_libfile))
    for lib_dir in default_library_dirs:
        try:
            ctypes.CDLL(os.path.join(lib_dir, target_libfile))
            target_libdir = lib_dir
        except OSError:
            pass
    return check_prefix(target_inc, target_libdir)


def check_for_dependencies(env, cfg, header, library):
    # First up: check in environment
    if env in os.environ:
        return get_location_from_env(env)
    # Next up, we try config file
    elif os.path.exists(cfg):
        return get_location_from_cfg(cfg)
    # Now we see if ctypes can help us
    if os.name == 'posix' or os.name == 'nt':
        target_inc, target_lib = get_location_from_ctypes(header, library)
    if None not in (target_inc, target_lib):
        print(
            "%s_LOCATION: %s found via ctypes in: %s, %s"
            % (env.split('_')[0], env.split('_')[0], target_inc, target_lib)
        )
        return (target_inc, target_lib)

    print("Reading %s location from %s failed." % (env.split('_')[0], cfg))
    print("Please place the base directory of your")
    print("%s install in %s and restart." % (env.split('_')[0], cfg))
    print("(ex: \"echo '/usr/local/' > %s\" )" % cfg)
    print("You can locate the path by looking for %s" % header)
    sys.exit(1)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('utilities', parent_package, top_path)
    config.add_subpackage("amr_kdtree")
    config.add_subpackage("poster")
    config.add_subpackage("answer_testing")
    config.add_subpackage("spatial")
    config.add_subpackage("grid_data_format")
    config.add_subpackage("parallel_tools")
    config.add_subpackage("lib")
    config.add_extension("data_point_utilities",
                         "yt/utilities/data_point_utilities.c",
                         libraries=["m"])
    config.add_subpackage("tests")
    config.add_subpackage("pyparselibconfig")
    config.make_config_py()  # installs __config__.py
    # config.make_svn_version_py()
    return config
