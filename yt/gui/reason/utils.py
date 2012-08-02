"""
Utilities for Reason

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from .bottle_mods import PayloadHandler
import base64
import types
import os

def load_script(filename):
    contents = open(filename).read()
    payload_handler = PayloadHandler()
    payload_handler.add_payload(
        {'type': 'script',
         'value': contents}
    )
    return

def deliver_image(im):
    if hasattr(im, 'read'):
        img_data = base64.b64encode(im.read())
    elif isinstance(im, types.StringTypes) and \
         im.endswith(".png"):
        img_data = base64.b64encode(open(im).read())
    elif isinstance(im, types.StringTypes):
        img_data = im
    else:
        raise RuntimeError
    ph = PayloadHandler()
    payload = {'type':'cell',
               'output': '',
               'input': '',
               'raw_input': '',
               'result_id': None,
               'image_data':img_data}
    ph.add_payload(payload)

def get_list_of_datasets():
    # Note that this instantiates the hierarchy.  This can be a costly
    # event.  However, we're going to assume that it's okay, if you have
    # decided to load up the parameter file.
    from yt.data_objects.static_output import _cached_pfs
    rv = []
    for fn, pf in sorted(_cached_pfs.items()):
        objs = []
        pf_varname = "_cached_pfs['%s']" % (fn)
        field_list = []
        if pf._instantiated_hierarchy is not None: 
            field_list = list(set(pf.h.field_list + pf.h.derived_field_list))
            field_list = [dict(text = f) for f in sorted(field_list)]
            for i,obj in enumerate(pf.h.objects):
                try:
                    name = str(obj)
                except ReferenceError:
                    continue
                objs.append(dict(name=name, type=obj._type_name,
                                 filename = '', field_list = [],
                                 varname = "%s.h.objects[%s]" % (pf_varname, i)))
        rv.append( dict(name = str(pf), children = objs, filename=fn,
                        type = "parameter_file",
                        varname = pf_varname, field_list = field_list) )
    return rv

def get_reasonjs_path():
    fn = "reason-js-20120623.zip"
    if "YT_DEST" not in os.environ:
        print
        print "*** You must set the environment variable YT_DEST ***"
        print "*** to point to the installation location!        ***"
        print
        raise IOError
    reasonjs_path = os.path.join(os.environ["YT_DEST"], "src", fn)
    if not os.path.isfile(reasonjs_path):
        print
        print "*** You are missing the Reason support files. You ***"
        print "*** You can get these by either rerunning the     ***"
        print "*** install script installing, or downloading     ***"
        print "*** them manually.                                ***"
        print "***                                               ***"
        print "*** FOR INSTANCE:                                 ***"
        print
        print "cd %s" % os.path.join(os.environ["YT_DEST"], "src")
        print "wget http://yt-project.org/dependencies/reason-js-20120623.zip"
        print
        raise IOError
    return reasonjs_path
