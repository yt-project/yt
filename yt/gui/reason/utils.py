"""
Utilities for Reason



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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
    # Note that this instantiates the index.  This can be a costly
    # event.  However, we're going to assume that it's okay, if you have
    # decided to load up the dataset.
    from yt.data_objects.static_output import _cached_datasets
    rv = []
    for fn, ds in sorted(_cached_datasets.items()):
        objs = []
        ds_varname = "_cached_datasets['%s']" % (fn)
        field_list = []
        if ds._instantiated_index is not None: 
            field_list = list(set(ds.field_list + ds.derived_field_list))
            field_list = [dict(text = f) for f in sorted(field_list)]
            for i,obj in enumerate(ds.h.objects):
                try:
                    name = str(obj)
                except ReferenceError:
                    continue
                objs.append(dict(name=name, type=obj._type_name,
                                 filename = '', field_list = [],
                                 varname = "%s.h.objects[%s]" % (ds_varname, i)))
        rv.append( dict(name = str(ds), children = objs, filename=fn,
                        type = "dataset",
                        varname = ds_varname, field_list = field_list) )
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
