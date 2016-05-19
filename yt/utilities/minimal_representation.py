"""
Skeleton objects that represent a few fundamental yt data types.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import abc
import json
import sys
from yt.utilities.on_demand_imports import _h5py as h5
import os
from uuid import uuid4
from yt.extern.six.moves import urllib
from yt.extern.six.moves import cPickle as pickle
from tempfile import TemporaryFile
from yt.config import ytcfg
from yt.funcs import \
    iterable, get_pbar, compare_dicts
from yt.extern.six import add_metaclass, string_types, b
from yt.utilities.exceptions import \
    YTHubRegisterError
from yt.utilities.logger import ytLogger as mylog
from yt.units.yt_array import \
    YTArray, \
    YTQuantity

if sys.version_info < (3, 0):
    from .poster.streaminghttp import register_openers
    from .poster.encode import multipart_encode
    register_openers()
else:
    # We don't yet have a solution for this, but it won't show up very often
    # anyway.
    pass


def _sanitize_list(flist):
    temp = []
    for item in flist:
        if isinstance(item, string_types):
            temp.append(b(item))
        elif isinstance(item, tuple) and \
                all(isinstance(i, string_types) for i in item):
            temp.append(tuple(_sanitize_list(list(item))))
        else:
            temp.append(item)
    return temp


def _serialize_to_h5(g, cdict):
    for item in cdict:
        if isinstance(cdict[item], (YTQuantity, YTArray)):
            g[item] = cdict[item].d
            g[item].attrs["units"] = str(cdict[item].units)
        elif isinstance(cdict[item], dict):
            _serialize_to_h5(g.create_group(item), cdict[item])
        elif cdict[item] is None:
            g[item] = "None"
        elif isinstance(cdict[item], list):
            g[item] = _sanitize_list(cdict[item])
        elif isinstance(cdict[item], tuple) and \
                all(isinstance(i, string_types) for i in cdict[item]):
            g[item] = tuple(_sanitize_list(cdict[item]))
        else:
            g[item] = cdict[item]


def _deserialize_from_h5(g, ds):
    result = {}
    for item in g:
        if item == "chunks":
            continue
        if "units" in g[item].attrs:
            if iterable(g[item]):
                result[item] = ds.arr(g[item][:], g[item].attrs["units"])
            else:
                result[item] = ds.quan(g[item][()],
                                       g[item].attrs["units"])
        elif isinstance(g[item], h5.Group):
            result[item] = _deserialize_from_h5(g[item], ds)
        elif g[item] == "None":
            result[item] = None
        else:
            try:
                result[item] = g[item][:]   # try array
            except ValueError:
                result[item] = g[item][()]  # fallback to scalar
    return result


class UploaderBar(object):
    pbar = None

    def __init__(self, my_name=""):
        self.my_name = my_name

    def __call__(self, name, prog, total):
        if self.pbar is None:
            self.pbar = get_pbar("Uploading %s " % self.my_name, total)
        self.pbar.update(prog)
        if prog == total:
            self.pbar.finish()


class ContainerClass(object):
    pass


@add_metaclass(abc.ABCMeta)
class MinimalRepresentation(object):

    def _update_attrs(self, obj, attr_list):
        for attr in attr_list:
            setattr(self, attr, getattr(obj, attr, None))
        if hasattr(obj, "ds"):
            self.output_hash = obj.ds._hash()
            self._ds_mrep = obj.ds._mrep
        if hasattr(obj, "data_source"):
            self.data_source_hash = obj.data_source._hash

    def __init__(self, obj):
        self._update_attrs(obj, self._attr_list)

    @abc.abstractmethod
    def _generate_post(self):
        pass

    @abc.abstractproperty
    def _attr_list(self):
        pass

    def _return_filtered_object(self, attrs):
        new_attrs = tuple(attr for attr in self._attr_list
                          if attr not in attrs)
        new_class = type('Filtered%s' % self.__class__.__name__,
                         (FilteredRepresentation,),
                         {'_attr_list': new_attrs})
        return new_class(self)

    @property
    def _attrs(self):
        return dict(((attr, getattr(self, attr)) for attr in self._attr_list))

    @classmethod
    def _from_metadata(cls, metadata):
        cc = ContainerClass()
        for a, v in metadata.values():
            setattr(cc, a, v)
        return cls(cc)

    def store(self, storage):
        if hasattr(self, "_ds_mrep"):
            self._ds_mrep.store(storage)
        metadata, (final_name, chunks) = self._generate_post()
        metadata['obj_type'] = self.type
        with h5.File(storage) as h5f:
            dset = str(uuid4())[:8]
            h5f.create_group(dset)
            _serialize_to_h5(h5f[dset], metadata)
            if len(chunks) > 0:
                g = h5f[dset].create_group('chunks')
                g.attrs['final_name'] = final_name
                for fname, fdata in chunks:
                    if isinstance(fname, (tuple, list)):
                        fname = "*".join(fname)

                    if isinstance(fdata, (YTQuantity, YTArray)):
                        g.create_dataset(fname, data=fdata.d,
                                         compression="lzf")
                        g[fname].attrs["units"] = str(fdata.units)
                    else:
                        g.create_dataset(fname, data=fdata, compression="lzf")

    def restore(self, storage, ds):
        pass

    def upload(self):
        api_key = ytcfg.get("yt", "hub_api_key")
        url = ytcfg.get("yt", "hub_url")
        if api_key == '':
            raise YTHubRegisterError
        metadata, (final_name, chunks) = self._generate_post()
        if hasattr(self, "_ds_mrep"):
            self._ds_mrep.upload()
        for i in metadata:
            if isinstance(metadata[i], np.ndarray):
                metadata[i] = metadata[i].tolist()
            elif hasattr(metadata[i], 'dtype'):
                metadata[i] = np.asscalar(metadata[i])
        metadata['obj_type'] = self.type
        if len(chunks) == 0:
            chunk_info = {'chunks': []}
        else:
            chunk_info = {'final_name': final_name, 'chunks': []}
            for cn, cv in chunks:
                chunk_info['chunks'].append((cn, cv.size * cv.itemsize))
        metadata = json.dumps(metadata)
        chunk_info = json.dumps(chunk_info)
        datagen, headers = multipart_encode({'metadata': metadata,
                                             'chunk_info': chunk_info,
                                             'api_key': api_key})
        request = urllib.request.Request(url, datagen, headers)
        # Actually do the request, and get the response
        try:
            rv = urllib.request.urlopen(request).read()
        except urllib.error.HTTPError as ex:
            if ex.code == 401:
                mylog.error("You must create an API key before uploading.")
                mylog.error("https://data.yt-project.org/getting_started.html")
                return
            else:
                raise ex
        uploader_info = json.loads(rv)
        new_url = url + "/handler/%s" % uploader_info['handler_uuid']
        for i, (cn, cv) in enumerate(chunks):
            f = TemporaryFile()
            np.save(f, cv)
            f.seek(0)
            pbar = UploaderBar("%s, % 2i/% 2i" %
                               (self.type, i + 1, len(chunks)))
            datagen, headers = multipart_encode({'chunk_data': f}, cb=pbar)
            request = urllib.request.Request(new_url, datagen, headers)
            rv = urllib.request.urlopen(request).read()

        datagen, headers = multipart_encode({'status': 'FINAL'})
        request = urllib.request.Request(new_url, datagen, headers)
        rv = json.loads(urllib.request.urlopen(request).read())
        mylog.info("Upload succeeded!  View here: %s", rv['url'])
        return rv

    def load(self, storage):
        return pickle.load(open(storage, 'r'))

    def dump(self, storage):
        with open(storage, 'w') as fh:
            pickle.dump(self, fh)


class FilteredRepresentation(MinimalRepresentation):

    def _generate_post(self):
        raise RuntimeError


class MinimalDataset(MinimalRepresentation):
    _attr_list = ("dimensionality", "refine_by", "domain_dimensions",
                  "current_time", "domain_left_edge", "domain_right_edge",
                  "unique_identifier", "current_redshift", "output_hash",
                  "cosmological_simulation", "omega_matter", "omega_lambda",
                  "hubble_constant", "name")
    type = 'simulation_output'

    def __init__(self, obj):
        super(MinimalDataset, self).__init__(obj)
        self.output_hash = obj._hash()
        self.name = str(obj)

    def _generate_post(self):
        metadata = self._attrs
        chunks = []
        return (metadata, (None, chunks))


class MinimalMappableData(MinimalRepresentation):

    _attr_list = ("field_data", "field", "weight_field", "axis", "output_hash",
                  "vm_type")

    def _generate_post(self):
        nobj = self._return_filtered_object(("field_data",))
        metadata = nobj._attrs
        chunks = [(arr, self.field_data[arr]) for arr in self.field_data]
        return (metadata, ('field_data', chunks))

    def _read_chunks(self, g, ds):
        for fname in g.keys():
            if '*' in fname:
                arr = tuple(fname.split('*'))
            else:
                arr = fname
            try:
                self.field_data[arr] = ds.arr(g[fname][:],
                                              g[fname].attrs["units"])
            except KeyError:
                self.field_data[arr] = g[fname][:]


class MinimalProjectionData(MinimalMappableData):
    type = 'proj'
    vm_type = "Projection"
    _attr_list = ("field_data", "field", "weight_field", "axis", "output_hash",
                  "center", "method", "field_parameters",
                  "data_source_hash")

    def restore(self, storage, ds):
        if hasattr(self, "_ds_mrep"):
            self._ds_mrep.restore(storage, ds)
        metadata, (final_name, chunks) = self._generate_post()
        with h5.File(storage, 'r') as h5f:
            for dset in h5f:
                stored_metadata = _deserialize_from_h5(h5f[dset], ds)
                if compare_dicts(metadata, stored_metadata):
                    self._read_chunks(h5f[dset]["chunks"], ds)
                    return True
        return False


class MinimalSliceData(MinimalMappableData):
    type = 'slice'
    vm_type = "Slice"
    weight_field = "None"


class MinimalImageCollectionData(MinimalRepresentation):
    type = "image_collection"
    _attr_list = ("name", "output_hash", "images", "image_metadata")

    def _generate_post(self):
        nobj = self._return_filtered_object(("images",))
        metadata = nobj._attrs
        chunks = [(fn, d) for fn, d in self.images]
        return (metadata, ('images', chunks))

_hub_categories = ("News", "Documents", "Simulation Management",
                   "Data Management", "Analysis and Visualization",
                   "Paper Repositories", "Astrophysical Utilities",
                   "yt Scripts")


class MinimalProjectDescription(MinimalRepresentation):
    type = "project"
    _attr_list = ("title", "url", "description", "category", "image_url")

    def __init__(self, title, url, description,
                 category, image_url=""):
        assert(category in _hub_categories)
        self.title = title
        self.url = url
        self.description = description
        self.category = category
        self.image_url = image_url

    def _generate_post(self):
        metadata = self._attrs
        return (metadata, ("chunks", []))


class MinimalNotebook(MinimalRepresentation):
    type = "notebook"
    _attr_list = ("title",)

    def __init__(self, filename, title=None):
        # First we read in the data
        if not os.path.isfile(filename):
            raise IOError(filename)
        self.data = open(filename).read()
        if title is None:
            title = json.loads(self.data)['metadata']['name']
        self.title = title
        self.data = np.fromstring(self.data, dtype='c')

    def _generate_post(self):
        metadata = self._attrs
        chunks = [("notebook", self.data)]
        return (metadata, ("chunks", chunks))


class ImageCollection(object):

    def __init__(self, ds, name):
        self.ds = ds
        self.name = name
        self.images = []
        self.image_metadata = []

    def add_image(self, fn, descr):
        self.image_metadata.append(descr)
        self.images.append((os.path.basename(fn), np.fromfile(fn, dtype='c')))
