import abc
import json
import os
from typing import Tuple
from uuid import uuid4

import numpy as np

from yt.funcs import compare_dicts, is_sequence
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.on_demand_imports import _h5py as h5py


def _sanitize_list(flist):
    temp = []
    for item in flist:
        if isinstance(item, str):
            temp.append(item.encode("latin-1"))
        elif isinstance(item, tuple) and all(isinstance(i, str) for i in item):
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
        elif isinstance(cdict[item], tuple) and all(
            isinstance(i, str) for i in cdict[item]
        ):
            g[item] = tuple(_sanitize_list(cdict[item]))
        else:
            g[item] = cdict[item]


def _deserialize_from_h5(g, ds):
    result = {}
    for item in g:
        if item == "chunks":
            continue
        if "units" in g[item].attrs:
            if is_sequence(g[item]):
                result[item] = ds.arr(g[item][:], g[item].attrs["units"])
            else:
                result[item] = ds.quan(g[item][()], g[item].attrs["units"])
        elif isinstance(g[item], h5py.Group):
            result[item] = _deserialize_from_h5(g[item], ds)
        elif g[item] == "None":
            result[item] = None
        else:
            try:
                result[item] = g[item][:]  # try array
            except ValueError:
                result[item] = g[item][()]  # fallback to scalar
    return result


class ContainerClass:
    pass


class MinimalRepresentation(metaclass=abc.ABCMeta):
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
        new_attrs = tuple(attr for attr in self._attr_list if attr not in attrs)
        new_class = type(
            f"Filtered{self.__class__.__name__}",
            (FilteredRepresentation,),
            {"_attr_list": new_attrs},
        )
        return new_class(self)

    @property
    def _attrs(self):
        return {attr: getattr(self, attr) for attr in self._attr_list}

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
        metadata["obj_type"] = self.type
        with h5py.File(storage, mode="r") as h5f:
            dset = str(uuid4())[:8]
            h5f.create_group(dset)
            _serialize_to_h5(h5f[dset], metadata)
            if len(chunks) > 0:
                g = h5f[dset].create_group("chunks")
                g.attrs["final_name"] = final_name
                for fname, fdata in chunks:
                    if isinstance(fname, (tuple, list)):
                        fname = "*".join(fname)

                    if isinstance(fdata, (YTQuantity, YTArray)):
                        g.create_dataset(fname, data=fdata.d, compression="lzf")
                        g[fname].attrs["units"] = str(fdata.units)
                    else:
                        g.create_dataset(fname, data=fdata, compression="lzf")

    def restore(self, storage, ds):
        pass

    def upload(self):
        raise NotImplementedError("This method hasn't been ported to python 3")

    def load(self, storage):
        raise NotImplementedError("This method hasn't been ported to python 3")

    def dump(self, storage):
        raise NotImplementedError("This method hasn't been ported to python 3")


class FilteredRepresentation(MinimalRepresentation):
    def _generate_post(self):
        raise RuntimeError


class MinimalDataset(MinimalRepresentation):
    _attr_list = (
        "dimensionality",
        "refine_by",
        "domain_dimensions",
        "current_time",
        "domain_left_edge",
        "domain_right_edge",
        "unique_identifier",
        "current_redshift",
        "output_hash",
        "cosmological_simulation",
        "omega_matter",
        "omega_lambda",
        "hubble_constant",
        "name",
    )
    type = "simulation_output"

    def __init__(self, obj):
        super().__init__(obj)
        self.output_hash = obj._hash()
        self.name = str(obj)

    def _generate_post(self):
        metadata = self._attrs
        chunks = []
        return (metadata, (None, chunks))


class MinimalMappableData(MinimalRepresentation):

    _attr_list: Tuple[str, ...] = (
        "field_data",
        "field",
        "weight_field",
        "axis",
        "output_hash",
        "vm_type",
    )

    def _generate_post(self):
        nobj = self._return_filtered_object(("field_data",))
        metadata = nobj._attrs
        chunks = [(arr, self.field_data[arr]) for arr in self.field_data]
        return (metadata, ("field_data", chunks))

    def _read_chunks(self, g, ds):
        for fname in g.keys():
            if "*" in fname:
                arr = tuple(fname.split("*"))
            else:
                arr = fname
            try:
                self.field_data[arr] = ds.arr(g[fname][:], g[fname].attrs["units"])
            except KeyError:
                self.field_data[arr] = g[fname][:]


class MinimalProjectionData(MinimalMappableData):
    type = "proj"
    vm_type = "Projection"
    _attr_list = (
        "field_data",
        "field",
        "weight_field",
        "axis",
        "output_hash",
        "center",
        "method",
        "field_parameters",
        "data_source_hash",
    )

    def restore(self, storage, ds):
        if hasattr(self, "_ds_mrep"):
            self._ds_mrep.restore(storage, ds)
        metadata, (final_name, chunks) = self._generate_post()
        with h5py.File(storage, mode="r") as h5f:
            for dset in h5f:
                stored_metadata = _deserialize_from_h5(h5f[dset], ds)
                if compare_dicts(metadata, stored_metadata):
                    self._read_chunks(h5f[dset]["chunks"], ds)
                    return True
        return False


class MinimalSliceData(MinimalMappableData):
    type = "slice"
    vm_type = "Slice"
    weight_field = "None"


class MinimalImageCollectionData(MinimalRepresentation):
    type = "image_collection"
    _attr_list = ("name", "output_hash", "images", "image_metadata")

    def _generate_post(self):
        nobj = self._return_filtered_object(("images",))
        metadata = nobj._attrs
        chunks = [(fn, d) for fn, d in self.images]
        return (metadata, ("images", chunks))


_hub_categories = (
    "News",
    "Documents",
    "Simulation Management",
    "Data Management",
    "Analysis and Visualization",
    "Paper Repositories",
    "Astrophysical Utilities",
    "yt Scripts",
)


class MinimalProjectDescription(MinimalRepresentation):
    type = "project"
    _attr_list = ("title", "url", "description", "category", "image_url")

    def __init__(self, title, url, description, category, image_url=""):
        assert category in _hub_categories
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
            raise OSError(filename)
        self.data = open(filename).read()
        if title is None:
            title = json.loads(self.data)["metadata"]["name"]
        self.title = title
        self.data = np.fromstring(self.data, dtype="c")

    def _generate_post(self):
        metadata = self._attrs
        chunks = [("notebook", self.data)]
        return (metadata, ("chunks", chunks))


class ImageCollection:
    def __init__(self, ds, name):
        self.ds = ds
        self.name = name
        self.images = []
        self.image_metadata = []

    def add_image(self, fn, descr):
        self.image_metadata.append(descr)
        self.images.append((os.path.basename(fn), np.fromfile(fn, dtype="c")))
