"""
YT is a package written primarily in python designed to make the task of
running Enzo  easier.  It contains facilities for creating Enzo data (currently
in prototype form) as well as runnning Enzo simulations, simulating the actions
of Enzo on various existing data, and analyzing output from Enzo in a
wide-variety of methods.

The changelog is available in
U{fancy<http://www.slac.stanford.edu/~mturk/yt_doc/changelog.html>} and
U{detailed<http://www.slac.stanford.edu/~mturk/yt_doc/changelog.txt>} formats.
These are automatically regenerated nightly.

An ever-growing HOWTO is also available at
U{http://www.slac.stanford.edu/~mturk/yt_doc/yt_howto.pdf} .

YT is divided into several packages, all named after characters from Snow
Crash.

Enki
====

    Enki is the package used to create data, and instantiate runs.  It supports
    creating Enzo Problems, and then using SWIG-interfaced Enzo calls to actually
    create the data for those problems.  Additionally, facilities are being
    developed to use Enki to directly execute runs.

Fido
====

    Fido is the messenger/protector of data.  It takes data outputs, puts them
    wherever you want, and then calls a function handler to deal with that data.
    Ultimately Fido will deal with all file-handling; submission of runs to a
    central (user-specific) database is in the works, and Fido will be the entity
    that moves the files in and out of storage.

Lagos
=====

    Lagos deals with data structures.  It defines things like EnzoGrid,
    EnzoData, Enzo2DData, EnzoSphere, etc.  If you want to handle actual data, use
    Lagos.

Raven
=====

    Raven is the interface to HippoDraw.  All data plotting goes through
    Raven.  An interface to Matplotlib is being worked on, and almost ready.
    (This includes interpolated plots as well as vector plots.  Also,
    contours.  And overlays.)

Deliverator
===========

    The Deliverator is a TurboGears-based system for querying and displaying
    images.  Images are dumped from Raven into local, web-accessible storage space,
    and then metadata about those images is submitted to The Deliverator.  The user
    (you) then goes to the Deliverator website and views those plots.

The base package YT provides facilities for configuration files and logging (via the
Python logger.)

Also, documentation should follow the epydoc format,
U{http://epydoc.sourceforge.net/} with the fields listed
U{here<http://epydoc.sourceforge.net/manual-fields.html>}.  I'm working on
making my way through the existing code to document it.  Symbols and other
markup information available
U{here<http://epydoc.sourceforge.net/api/epydoc.markup.epytext-module.html>}.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@see: U{Snow Crash<http://en.wikipedia.org/wiki/Snow_Crash>}
@see: U{Enki<http://en.wikipedia.org/wiki/Enki>}
"""

import importer
importer.install()

from logger import *
from config import *
from funcs import *
