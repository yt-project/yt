# -*- coding: utf-8 -*-
"""
    yt.extensions
    ~~~~~~~~~~~~~

    Redirect imports for extensions.  This module basically makes it possible
    for us to transition from ytext.foo to ytext_foo without having to
    force all extensions to upgrade at the same time.

    When a user does ``from yt.extensions.foo import bar`` it will attempt to
    import ``from yt_foo import bar`` first and when that fails it will
    try to import ``from ytext.foo import bar``.

    We're switching from namespace packages because it was just too painful for
    everybody involved.

    :copyright: (c) 2015 by Armin Ronacher.
    :license: BSD, see LICENSE for more details.
"""

# This source code is originally from flask, in the flask/ext/__init__.py file.


def setup():
    from ..exthook import ExtensionImporter
    importer = ExtensionImporter(['yt_%s', 'ytext.%s'], __name__)
    importer.install()


setup()
del setup
